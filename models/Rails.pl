% UAV properties
1.0::standard; 0.0::special.
initial_charge ~ normal(90, 5).
charge_cost ~ normal(-0.1, 0.2).
weight ~ normal(2.0, 0.1).

1/10::fog; 9/10::clear.

% Sufficient charge is defined to be
% enough to return to the operator
can_return(R, C) :-
    B is initial_charge,
    O is charge_cost,
    D is distance(R, C, operator),
    0 < B + (2 * O * D).

% Visual line of sight
vlos(R, C) :- 
    fog, distance(R, C, operator) < 250;
    clear, distance(R, C, operator) < 500.

% Simplified OPEN flight category
open_flight(R, C) :- 
    standard, vlos(R, C), weight < 25.


fly_eastward_of_primary(R, C) :- 
    distance_east(R, C, primary) < 100, 
    distance_east(R, C, primary) > 5.

% permit to fly eastward of primary but not between primary road
permit(R, C) :-
    fly_eastward_of_primary(R, C),
    \+ betweenCustom(R, C, primary).

% Stay close to rails and keep enough battery charge to return
landscape(R, C) :-  can_return(R, C), open_flight(R, C), permit(R, C).
