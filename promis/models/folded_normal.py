"""This module includes an abstraction of folded normal distributions."""

#
# Copyright (c) Simon Kohaut, Honda Research Institute Europe GmbH
#
# This file is part of ProMis and licensed under the BSD 3-Clause License.
# You should have received a copy of the BSD 3-Clause License along with ProMis.
# If not, see https://opensource.org/license/bsd-3-clause/.
#

# Standard Library
from typing import cast

# Third Party
from numpy import ndarray
from scipy.stats import foldnorm


class FoldedNormal:

    """A folded normal distribution.

    Args:
        mean: The mean of the distribution 
        covariance: The covariance of the distribution

    References:
        - https://en.wikipedia.org/wiki/Folded_normal_distribution
    """

    def __init__(self, mean: float, variance: float):
        # Sanity checks on given parameters
        assert variance > 0 # variance has to be positive

        # Assign values
        self.mean = mean
        self.variance = variance
        std_div = variance**0.5

        # Abstract away from the scipy implementation
        self.distribution = foldnorm(mean/std_div, 0, std_div)

    @property
    def x(self) -> ndarray:
        return self.mean

    @property
    def P(self) -> ndarray:
        return self.variance

    def sample(self, number_of_samples: int = 1) -> ndarray:
        """Draw a number of samples following this folded normal's distribution.

        Args:
            number_of_samples: The number of samples to draw

        Returns:
            The drawn samples
        """

        assert number_of_samples >= 1, "Number of samples cannot be negative or zero!"

        
        samples = self.distribution.rvs(size=number_of_samples)

        return cast(ndarray, samples)

    def cdf(self, x: ndarray) -> float:
        """Compute the CDF as integral of the PDF from negative infinity up to x.

        Args:
            x: The upper bound of the integral

        Returns:
            The probability of a value being less than x
        """

        return cast(float, self.distribution.cdf(x))

    def __call__(self, value: float) -> float:
        """Evaluate the folded normal at the given location, i.e. obtain the probability density.

        Args:
            value: Where to evaluate the folded normal, of dimension ``(n, 1)``

        Returns:
            The probability density at the given location
        """

        # Compute and return weighted probability density function
        return cast(float, self.distribution.pdf(value))

    def __eq__(self, other) -> bool:
        """Checks if two folded normal distributions are equal.

        Args:
            other: The distribution to compare within

        Returns:
            Whether the two distributions are the same
        """

        return (
            cast(bool, (self.mean == other.mean))
            and cast(bool, (self.covariance == other.covariance))
        )
