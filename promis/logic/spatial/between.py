"""This module implements a distributional predicate of distances to sets of map features."""

#
# Copyright (c) Simon Kohaut, Honda Research Institute Europe GmbH
#
# This file is part of ProMis and licensed under the BSD 3-Clause License.
# You should have received a copy of the BSD 3-Clause License along with ProMis.
# If not, see https://opensource.org/license/bsd-3-clause/.
#

# Standard Library
from itertools import product
from pathlib import Path
from typing import cast
import math

# Third Party
from numpy import mean, unravel_index
from shapely.strtree import STRtree
from shapely import shortest_line
from shapely.geometry import LineString, Point, Polygon
# ProMis
from promis.geo import CartesianLocation, CartesianMap, LocationType, PolarMap, RasterBand


class Between:

    """TODO"""

    def __init__(
        self,
        probability: RasterBand,
        location_type: LocationType,
    ) -> None:
        # Setup attributes
        self.probability = probability
        self.location_type = location_type

    def save_as_plp(self, path: Path) -> None:
        with open(path, "w") as plp_file:
            plp_file.write(self.to_distributional_clauses())

    def to_distributional_clauses(self) -> str:
        code = ""
        for index in product(
            range(self.probability.data.shape[0]), range(self.probability.data.shape[1])
        ):
            code += self.index_to_distributional_clause(index)

        return code

    def index_to_distributional_clause(self, index: tuple[int, int]) -> str:
        feature_name = self.location_type.name.lower()

        relation = f"between(row_{index[1]}, column_{index[0]}, {feature_name}).\n"

        if self.probability.data[index] == 1.0:
            return relation
        else:
            return f"{self.probability.data[index]}::{relation}"

    def split(self) -> "list[list[Between]] | Between":
        probability_splits = self.probability.split()

        if isinstance(probability_splits, RasterBand):
            return self

        return [
            [
                Between(probability_splits[0][0], self.location_type),
                Between(probability_splits[0][1], self.location_type),
            ],
            [
                Between(probability_splits[1][0], self.location_type),
                Between(probability_splits[1][1], self.location_type),
            ],
        ]

    @classmethod
    def from_map(
        cls,
        map_: PolarMap | CartesianMap,
        location_type: LocationType,
        resolution: tuple[int, int],
        number_of_samples: int = 50,
        radius: int = 300,
    ) -> "Between":
        # Setup attributes
        cartesian_map = map_ if isinstance(map_, CartesianMap) else map_.to_cartesian()

        # If map is empy return
        if cartesian_map.features is None:
            return None

        # Get all relevant features
        features = [
            feature for feature in cartesian_map.features if feature.location_type == location_type
        ]
        if not features:
            return None

        # Construct an STR tree per collection of varitions of features
        str_trees = [
            STRtree(
                [
                    feature.sample()[0].geometry
                    if feature.distribution is not None
                    else feature.geometry
                    for feature in features
                ]
            )
            for _ in range(number_of_samples)
        ]

        # Prepare raster bands
        probability = RasterBand(
            resolution, cartesian_map.origin, cartesian_map.width, cartesian_map.height
        )

        # Compute parameters of normal distributions for each location
        for i, location in enumerate(probability.cartesian_locations.values()):
            index = unravel_index(i, probability.data.shape)
            probability.data[index] = cls.compute_probabilities(location, str_trees, radius)

        # Create and return Over object
        return cls(probability, location_type)

    @staticmethod
    def compute_probabilities(location: CartesianLocation, str_trees: list[STRtree], radius: int) -> float:
        """Computes the probability for a location to be between geometries of some type.

        Args:
            location: The location to compute the  probability for
            str_trees: Random variations of the features of a map indexible by an STRtree each

        Returns:
            The probability for a location to be between geometries of some type
        """
        probabilities = []
        for str_tree in str_trees:
            indices = str_tree.query(location.geometry, 'dwithin', radius)
            features = str_tree.geometries.take(indices).tolist()
            probability = Between.compute_probability(location, features)
            probabilities.append(probability)

        return cast(float, mean(probabilities))
            

    @staticmethod
    def compute_probability(location: CartesianLocation, features: list[Point | LineString | Polygon]) -> float:
        """Computes the value for a location to be between geometries of some type.

        Args:
            location: The location to compute the  probability for
            features: list of feature concerning the computation of the same type

        Returns:
            The value for a location to be between geometries of some type
        """
        
        distances_E = []
        distances_N = []
        for feature in features:
            line = shortest_line(
                location.geometry,
                feature
            )
            distances_E.append(
                line.coords[1][0] - line.coords[0][0]         # des_x - location_x
            )
            distances_N.append(
                line.coords[1][1] - line.coords[0][1]         # des_y - location_y
            )

        num_of_pos_neg_E = (sum(map(lambda x : x > 0, distances_E)), sum(map(lambda x : x < 0, distances_E)))
        num_of_pos_neg_N = (sum(map(lambda x : x > 0, distances_N)), sum(map(lambda x : x < 0, distances_N)))

        prob_E = 0
        prob_N = 0

        def check_zero(pos_neg):
            return pos_neg[0] != 0 and pos_neg[1] != 0
        
        if check_zero(num_of_pos_neg_E):
            prob_E = num_of_pos_neg_E[1] / num_of_pos_neg_E[0] if num_of_pos_neg_E[0] > num_of_pos_neg_E[1] else num_of_pos_neg_E[0] / num_of_pos_neg_E[1]
        if check_zero(num_of_pos_neg_N):
            prob_N = num_of_pos_neg_N[1] / num_of_pos_neg_N[0] if num_of_pos_neg_N[0] > num_of_pos_neg_N[1] else num_of_pos_neg_N[0] / num_of_pos_neg_N[1]
        
        return cast(float, prob_E * prob_N)
    
    @staticmethod
    def compute_probability_alter(location: CartesianLocation, features: list[Point | LineString | Polygon]) -> float:
        """Computes the value for a location to be between geometries of some type. alternative
        where we compare 2 vector of the nearest feature with its angle and magnitude.

        Args:
            location: The location to compute the  probability for
            features: list of feature concerning the computation of the same type

        Returns:
            The value for a location to be between geometries of some type
        """
        
        # create str_tree
        str_tree = STRtree(features)
        # get the shortest line to the nearest obj
        nearest_geo = str_tree.geometries.take(str_tree.nearest(location.geometry))
        first_nearest_line = shortest_line(
            location.geometry,
            nearest_geo
        )
        nearest_vector = [first_nearest_line.coords[1][0] - first_nearest_line.coords[0][0],
                          first_nearest_line.coords[1][1] - first_nearest_line.coords[0][1]]
        # remove the nearest geo
        features_remove = [feature for feature in features if feature is not nearest_geo]
        str_tree = STRtree(features_remove)
        second_nearest_line = shortest_line(
            location.geometry,
            str_tree.geometries.take(str_tree.nearest(location.geometry))
        )
        sec_nearest_vector = [second_nearest_line.coords[1][0] - second_nearest_line.coords[0][0],
                              second_nearest_line.coords[1][1] - second_nearest_line.coords[0][1]]
        def angle_between_vectors(u, v):
            dot_product = sum(i*j for i, j in zip(u, v))
            norm_u = math.sqrt(sum(i**2 for i in u))
            norm_v = math.sqrt(sum(i**2 for i in v))
            cos_theta = dot_product / (norm_u * norm_v)
            angle_rad = math.acos(cos_theta)
            angle_deg = math.degrees(angle_rad)
            return angle_deg, norm_u, norm_v
        # calculate angle between 2 vector
        angle, length_first_line, length_second_line = angle_between_vectors(nearest_vector, sec_nearest_vector)
        if (angle < 120 or length_first_line == 0 or length_second_line == 0):
            return 0
        
        result = length_first_line / length_second_line if length_first_line > length_second_line else length_second_line / length_first_line

        return result