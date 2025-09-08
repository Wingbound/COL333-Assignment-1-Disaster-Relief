#include "solver.h"
#include <bits/stdc++.h>
#include <chrono>
#include <unordered_map>
#include <queue>
#include <unordered_set>
#include <random>

using namespace std;
// You can add any helper functions or classes you need here.

struct VillageAlt{
    Village village;
    int dry_food_recieved;
    int perishable_food_recieved;
    int other_supplies_recieved;
    bool is_fulfilled;
    int food_required() const {
        return village.population * 9 - (dry_food_recieved + perishable_food_recieved);
    }
    int other_required() const {
        return village.population - other_supplies_recieved;
    }
};

struct HelicopterAlt{
    Helicopter helicopter;
    double cumulative_distance_travelled;
    double trip_distance;
};

struct PQByHeuristic {
    bool operator()(const pair<VillageAlt*, double>& a,
                    const pair<VillageAlt*, double>& b) const {
        return a.second < b.second; // max-heap based on heuristic value
    }
};

bool satisfaction(vector<VillageAlt>& villages_alt){
    bool all_good = true;
    for(VillageAlt& v: villages_alt){
        all_good = all_good && v.is_fulfilled;
    }
    return all_good;
}

bool time_check(chrono::high_resolution_clock::time_point start, int time_limit_minutes) {
    auto now = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> elapsed = now - start;
    return elapsed.count() < time_limit_minutes * 60 * 1000 - 35*1000;
}

// priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, comparator> unvisited_villages_queue(vector<VillageAlt>& villages_alt, Point curr, unordered_set<int>& visited_villages_ids){
//     priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, comparator> pq;

//     for(VillageAlt& v : villages_alt){
//         if(visited_villages_ids.count(v.village.id)) continue;

//         double dist = distance(curr, v.village.coords);
//         pq.push(make_pair(&v, dist));
//     }

//     return pq;
// }

string vecToString(vector<HelicopterAlt>& helicopters_permutations) {
    string s;
    for (HelicopterAlt h : helicopters_permutations) {
        s += to_string(h.helicopter.id);
    }
    return s;
}

double calculateTripDist(const Trip& trip, const Helicopter& heli, const ProblemData& problem) {
    
    double total_distance = 0.0;
    Point current_location = problem.cities[heli.home_city_id - 1 ]; 
    for (const Drop& drop : trip.drops) {
        Point village_coords = problem.villages[drop.village_id - 1].coords; 
        total_distance += distance(current_location, village_coords);
        current_location = village_coords;
    }

    // Return to home city
    total_distance += distance(current_location, problem.cities[heli.home_city_id - 1]); 

    return total_distance;
}

double calculateTripCost(const Helicopter& helicopter, const Trip& trip, const ProblemData& problem) {
    double trip_distance = calculateTripDist(trip, helicopter, problem);
    return helicopter.fixed_cost + helicopter.alpha * trip_distance;
}

// double calculateTripValue(const Trip& trip, const ProblemData& problem, const Helicopter& helicopter) {
//     double total_value = 0.0;
//     for (const Drop& drop : trip.drops) {
//         total_value += drop.dry_food * problem.packages[0].value;
//         total_value += drop.perishable_food * problem.packages[1].value;
//         total_value += drop.other_supplies * problem.packages[2].value;
//     }
//     double trip_cost = calculateTripCost(helicopter, trip, problem);
//     return total_value - trip_cost;
// }

// double evaluatePlan(const ProblemData& problem, const HelicopterPlan& hp) {
//     double total_value = 0.0;

//     for (const Trip& trip : hp.trips) {
//         const Helicopter& helicopter = problem.helicopters[hp.helicopter_id - 1];
//         total_value += calculateTripValue(trip, problem, helicopter);
//     }

//     return total_value;
// }

// double evaluateSolution(const ProblemData& problem, const Solution& solution) {
//     double total_value = 0.0;
//     double total_cost = 0.0;
//     for (const HelicopterPlan& hp: solution) {
//         const Helicopter& helicopter = problem.helicopters[hp.helicopter_id - 1]; 
//         for (const Trip& trip : hp.trips) {
//             total_cost += calculateTripCost(helicopter, trip, problem);
//         }
//     }


//     return total_value - total_cost;

// }

// Helper function to calculate total supplies delivered to each village across all helicopters
map<int, array<double, 3>> calculateVillageSupplies(const ProblemData& problem, const Solution& solution) {
    map<int, array<double, 3>> village_supplies; // village_id -> [dry_food, perishable_food, other_supplies]
    
    // Initialize all villages with zero supplies
    for (const Village& village : problem.villages) {
        village_supplies[village.id] = {0.0, 0.0, 0.0};
    }
    
    // Accumulate supplies from all helicopters and trips
    for (const HelicopterPlan& hp : solution) {
        for (const Trip& trip : hp.trips) {
            for (const Drop& drop : trip.drops) {
                village_supplies[drop.village_id][0] += drop.dry_food;
                village_supplies[drop.village_id][1] += drop.perishable_food;
                village_supplies[drop.village_id][2] += drop.other_supplies;
            }
        }
    }
    
    return village_supplies;
}

// Calculate value for a single village based on diminishing returns
double calculateVillageValue(const Village& village, const array<double, 3>& supplies_delivered, 
                           const ProblemData& problem) {
    double total_value = 0.0;
    
    // Food supplies (dry + perishable) - capped at 9 times number of people
    double total_food_delivered = supplies_delivered[0] + supplies_delivered[1];
    double max_useful_food = 9.0 * village.population;
    double useful_food = min(total_food_delivered, max_useful_food);
    
    // Calculate value for useful food (proportional split between dry and perishable)
    if (total_food_delivered > 0) {
        double dry_ratio = supplies_delivered[0] / total_food_delivered;
        double perishable_ratio = supplies_delivered[1] / total_food_delivered;
        
        total_value += useful_food * dry_ratio * problem.packages[0].value;
        total_value += useful_food * perishable_ratio * problem.packages[1].value;
    }
    
    // Other supplies - capped at 1 times number of people
    double max_useful_other = 1.0 * village.population;
    double useful_other = min(supplies_delivered[2], max_useful_other);
    total_value += useful_other * problem.packages[2].value;
    
    return total_value;
}

// Main evaluation function with diminishing returns
double evaluateSolution(const ProblemData& problem, const Solution& solution) {
    double total_value = 0.0;
    double total_cost = 0.0;
    
    // Calculate total costs from all helicopter operations
    for (const HelicopterPlan& hp : solution) {
        const Helicopter& helicopter = problem.helicopters[hp.helicopter_id - 1];
        for (const Trip& trip : hp.trips) {
            total_cost += calculateTripCost(helicopter, trip, problem);
        }
    }
    
    // Calculate total value with diminishing returns per village
    map<int, array<double, 3>> village_supplies = calculateVillageSupplies(problem, solution);
    
    for (const Village& village : problem.villages) {
        array<double, 3> supplies = village_supplies[village.id];
        total_value += calculateVillageValue(village, supplies, problem);
    }
    
    return total_value - total_cost;
}

// Alternative evaluation that provides more detailed breakdown (useful for debugging)
struct EvaluationBreakdown {
    double total_value;
    double total_cost;
    double net_score;
    map<int, double> village_values;
    map<int, array<double, 3>> village_supplies_delivered;
    map<int, array<double, 3>> village_supplies_wasted; // supplies beyond useful limits
};

EvaluationBreakdown evaluateSolutionDetailed(const ProblemData& problem, const Solution& solution) {
    EvaluationBreakdown breakdown;
    breakdown.total_value = 0.0;
    breakdown.total_cost = 0.0;
    
    // Calculate costs
    for (const HelicopterPlan& hp : solution) {
        const Helicopter& helicopter = problem.helicopters[hp.helicopter_id - 1];
        for (const Trip& trip : hp.trips) {
            breakdown.total_cost += calculateTripCost(helicopter, trip, problem);
        }
    }
    
    // Calculate supplies delivered to each village
    breakdown.village_supplies_delivered = calculateVillageSupplies(problem, solution);
    
    // Calculate value and waste for each village
    for (const Village& village : problem.villages) {
        array<double, 3> supplies = breakdown.village_supplies_delivered[village.id];
        
        // Calculate value
        double village_value = calculateVillageValue(village, supplies, problem);
        breakdown.village_values[village.id] = village_value;
        breakdown.total_value += village_value;
        
        // Calculate waste
        double total_food = supplies[0] + supplies[1];
        double max_useful_food = 9.0 * village.population;
        double food_waste = max(0.0, total_food - max_useful_food);

        double max_useful_other = 1.0 * village.population;
        double other_waste = max(0.0, supplies[2] - max_useful_other);
        
        // Distribute food waste proportionally
        if (total_food > 0 && food_waste > 0) {
            double dry_ratio = supplies[0] / total_food;
            double perishable_ratio = supplies[1] / total_food;
            breakdown.village_supplies_wasted[village.id] = {
                food_waste * dry_ratio,
                food_waste * perishable_ratio,
                other_waste
            };
        } else {
            breakdown.village_supplies_wasted[village.id] = {0.0, 0.0, other_waste};
        }
    }
    
    breakdown.net_score = breakdown.total_value - breakdown.total_cost;
    return breakdown;
}

// Evaluate a single helicopter plan (useful for local search)
double evaluateHelicopterPlan(const ProblemData& problem, const HelicopterPlan& hp, 
                             const Solution& full_solution) {
    // This function calculates the contribution of a single helicopter plan
    // considering the diminishing returns based on what other helicopters have already delivered
    
    // Get supplies from other helicopters
    map<int, array<double, 3>> other_supplies;
    for (const Village& village : problem.villages) {
        other_supplies[village.id] = {0.0, 0.0, 0.0};
    }
    
    for (const HelicopterPlan& other_hp : full_solution) {
        if (other_hp.helicopter_id == hp.helicopter_id) continue; // Skip the current helicopter
        
        for (const Trip& trip : other_hp.trips) {
            for (const Drop& drop : trip.drops) {
                other_supplies[drop.village_id][0] += drop.dry_food;
                other_supplies[drop.village_id][1] += drop.perishable_food;
                other_supplies[drop.village_id][2] += drop.other_supplies;
            }
        }
    }
    
    double plan_value = 0.0;
    double plan_cost = 0.0;
    
    // Calculate cost of this helicopter's operations
    const Helicopter& helicopter = problem.helicopters[hp.helicopter_id - 1];
    for (const Trip& trip : hp.trips) {
        plan_cost += calculateTripCost(helicopter, trip, problem);
    }
    
    // Calculate marginal value added by this helicopter
    for (const Trip& trip : hp.trips) {
        for (const Drop& drop : trip.drops) {
            // Current supplies at this village (from other helicopters)
            array<double, 3> current_supplies = other_supplies[drop.village_id];
            
            // Supplies after adding this drop
            array<double, 3> new_supplies = {
                current_supplies[0] + drop.dry_food,
                current_supplies[1] + drop.perishable_food,
                current_supplies[2] + drop.other_supplies
            };
            
            // Find the village
            const Village* village = nullptr;
            for (const Village& v : problem.villages) {
                if (v.id == drop.village_id) {
                    village = &v;
                    break;
                }
            }
            
            if (village) {
                // Calculate marginal value
                double value_before = calculateVillageValue(*village, current_supplies, problem);
                double value_after = calculateVillageValue(*village, new_supplies, problem);
                plan_value += (value_after - value_before);
                
                // Update other_supplies for next drops to same village
                other_supplies[drop.village_id] = new_supplies;
            }
        }
    }
    
    return plan_value - plan_cost;
}

// Utility function to print evaluation breakdown (useful for debugging)
void printEvaluationBreakdown(const EvaluationBreakdown& breakdown) {
    cout << "=== Evaluation Breakdown ===" << endl;
    cout << "Total Value: " << breakdown.total_value << endl;
    cout << "Total Cost: " << breakdown.total_cost << endl;
    cout << "Net Score: " << breakdown.net_score << endl;
    cout << "\nPer Village:" << endl;
    
    for (const auto& [village_id, value] : breakdown.village_values) {
        auto supplies = breakdown.village_supplies_delivered.at(village_id);
        auto waste = breakdown.village_supplies_wasted.at(village_id);
        
        cout << "Village " << village_id << ": Value=" << value 
             << " | Delivered: [" << supplies[0] << ", " << supplies[1] << ", " << supplies[2] << "]"
             << " | Wasted: [" << waste[0] << ", " << waste[1] << ", " << waste[2] << "]" << endl;
    }
}

bool isFeasibleTrip(const ProblemData& problem, const HelicopterPlan& hp, const Trip& trip) {
    // Check for weight and distance constraints
    double total_weight = 0.0;
    for (const Drop& drop : trip.drops) {
        total_weight += drop.dry_food * problem.packages[0].weight + 
        drop.perishable_food * problem.packages[1].weight + drop.other_supplies * problem.packages[2].weight;
    }
    if (total_weight > problem.helicopters[hp.helicopter_id - 1].weight_capacity) {
        return false;
    }

    double trip_distance = calculateTripDist(trip, problem.helicopters[hp.helicopter_id - 1], problem);
    if (trip_distance > problem.helicopters[hp.helicopter_id - 1].distance_capacity) {
        return false;
    }

    return true;
}

bool isFeasible(const ProblemData& problem, const HelicopterPlan& hp) {
    // check for weight and distance of trips and cumulative distance for a heli
    for (const Trip& trip : hp.trips) {
        if (!isFeasibleTrip(problem, hp, trip)) {
            cout << "trip is not feasible" << endl;
            return false;
        }
    }

    // Check cumulative distance for the helicopter
    double total_distance = 0.0;
    for (const Trip& trip : hp.trips) {
        total_distance += calculateTripDist(trip, problem.helicopters[hp.helicopter_id - 1], problem);
    }
    if (total_distance > problem.d_max) {
        cout << "cumulative distance exceeded for helicopter " << endl;
        return false;
    }

    return true;
}

HelicopterPlan getCutPasteNeighbor(const ProblemData& problem, const HelicopterPlan& hp) {
    int helicopter_id = hp.helicopter_id;
    vector<Trip> trips = hp.trips;
    if (trips.empty()) return hp; // No trips to modify
    random_device rd;
    mt19937 gen(rd());
    int num_trips = trips.size();
    if (num_trips < 2) {
        cout << "didnt change";
        return hp;
    }
    uniform_int_distribution<> dis(0, num_trips - 1);
    int trip1 = dis(gen);
    // trip2 chosen as trip with leat total distance
    double min_dist = 1e18;
    int trip2 = -1;
    for (int i = 0; i <= (int) trips.size() - 1; i++) {
        if (i == trip1) continue;
        double x = calculateTripDist(trips[i], problem.helicopters[helicopter_id - 1], problem);
        if (x < min_dist) {
            trip2 = i;
            min_dist = x;
        } 
    }
    // int trip2;
    // do {
    //     trip2 = dis(gen);
    // } while (trip2 == trip1 || trips[trip2].drops.empty() || trips.size() < 2); // Ensure trip2 is different and has drops
    // choosing the least valuable drop, according to the max distance travelled for that drop
    int least_valuable_drop = -1;
    double score = 0.0; 
    double least_score = -1e18;
    Point prev = problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 - 1 ];
    for (int i = 0; i < trips[trip1].drops.size(); i++) {
        int village_id = trips[trip1].drops[i].village_id;
        Point curr = problem.villages[village_id - 1].coords;
        score = distance(curr, prev);
        if (score > least_score) { // currently judging only on the basis of distance
            least_valuable_drop = i;
            least_score = score;
        }
        prev = curr;
    }
    if (trip2 == -1 || least_valuable_drop == -1 || trips[trip1].drops.empty()) {
        cout << "didnt change " << trip2 << " " << least_valuable_drop << " " << trips[trip1].drops.size() << endl;

        return hp;
    }
    // remove the least valuable drop from trip1 and plug into trip2

    HelicopterPlan new_hp = hp; // Deep copy
    Drop drop_to_move = new_hp.trips[trip1].drops[least_valuable_drop];
    new_hp.trips[trip1].drops.erase(new_hp.trips[trip1].drops.begin() + least_valuable_drop);

    // Smart insertion into trip2: find best position to minimize added distance
    Point drop_point = problem.villages[drop_to_move.village_id - 1].coords;
    int best_pos = 0;
    double min_added_dist = 1e18;
    const vector<Drop>& trip2_drops = new_hp.trips[trip2].drops;

    for (int i = 0; i <= trip2_drops.size(); i++) {
        Point prev = (i == 0)
            ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 - 1]
            : problem.villages[trip2_drops[i - 1].village_id - 1].coords;

        Point next = (i == trip2_drops.size())
            ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 - 1]
            : problem.villages[trip2_drops[i].village_id - 1].coords;

        double added_dist = distance(prev, drop_point) + distance(drop_point, next) - distance(prev, next);
        if (added_dist < min_added_dist) {
            min_added_dist = added_dist;
            best_pos = i;
        }
    }

    new_hp.trips[trip2].drops.insert(new_hp.trips[trip2].drops.begin() + best_pos, drop_to_move);
    new_hp.trips[trip2].dry_food_pickup += drop_to_move.dry_food;
    new_hp.trips[trip2].perishable_food_pickup += drop_to_move.perishable_food;
    new_hp.trips[trip2].other_supplies_pickup += drop_to_move.other_supplies;

    // final check to see if new_hp is feasible

    if (!isFeasible(problem, new_hp)) {
        cout << "New plan is infeasible, returning original plan.\n";
        return hp; // Return original plan if new one is infeasible
    }

    return new_hp;
}


HelicopterPlan getSwapNeighbor(const ProblemData& problem, const HelicopterPlan& hp) {
    int helicopter_id = hp.helicopter_id;
    vector<Trip> trips = hp.trips;
    
    if (trips.empty()) return hp;
    
    random_device rd;
    mt19937 gen(rd());
    
    // Find two trips that both have drops
    vector<int> valid_trips;
    for (int i = 0; i < trips.size(); i++) {
        if (!trips[i].drops.empty()) {
            valid_trips.push_back(i);
        }
    }
    
    if (valid_trips.size() < 2) {
        cout << "Not enough trips with drops for swap" << endl;
        return hp;
    }
    
    // Select two different trips
    uniform_int_distribution<> trip_dis(0, valid_trips.size() - 1);
    int trip1_idx = trip_dis(gen);
    int trip2_idx;
    do {
        trip2_idx = trip_dis(gen);
    } while (trip2_idx == trip1_idx);
    
    int trip1 = valid_trips[trip1_idx];
    int trip2 = valid_trips[trip2_idx];
    
    // Select random drops from each trip
    uniform_int_distribution<> drop1_dis(0, trips[trip1].drops.size() - 1);
    uniform_int_distribution<> drop2_dis(0, trips[trip2].drops.size() - 1);
    int drop1_idx = drop1_dis(gen);
    int drop2_idx = drop2_dis(gen);
    
    HelicopterPlan new_hp = hp;
    
    // Get the drops to swap
    Drop drop1 = new_hp.trips[trip1].drops[drop1_idx];
    Drop drop2 = new_hp.trips[trip2].drops[drop2_idx];
    
    // Remove both drops
    new_hp.trips[trip1].drops.erase(new_hp.trips[trip1].drops.begin() + drop1_idx);
    new_hp.trips[trip2].drops.erase(new_hp.trips[trip2].drops.begin() + drop2_idx);
    
    // Update trip1 supplies (remove drop1, will add drop2)
    new_hp.trips[trip1].dry_food_pickup -= drop1.dry_food;
    new_hp.trips[trip1].perishable_food_pickup -= drop1.perishable_food;
    new_hp.trips[trip1].other_supplies_pickup -= drop1.other_supplies;
    
    // Update trip2 supplies (remove drop2, will add drop1)
    new_hp.trips[trip2].dry_food_pickup -= drop2.dry_food;
    new_hp.trips[trip2].perishable_food_pickup -= drop2.perishable_food;
    new_hp.trips[trip2].other_supplies_pickup -= drop2.other_supplies;
    
    // Find best position for drop2 in trip1
    Point drop2_point = problem.villages[drop2.village_id].coords;
    int best_pos1 = 0;
    double min_added_dist1 = 1e18;
    
    for (int i = 0; i <= new_hp.trips[trip1].drops.size(); i++) {
        Point prev = (i == 0)
            ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1]
            : problem.villages[new_hp.trips[trip1].drops[i - 1].village_id].coords;
            
        Point next = (i == new_hp.trips[trip1].drops.size())
            ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1]
            : problem.villages[new_hp.trips[trip1].drops[i].village_id].coords;
            
        double added_dist = distance(prev, drop2_point) + distance(drop2_point, next) - distance(prev, next);
        if (added_dist < min_added_dist1) {
            min_added_dist1 = added_dist;
            best_pos1 = i;
        }
    }
    
    // Find best position for drop1 in trip2
    Point drop1_point = problem.villages[drop1.village_id].coords;
    int best_pos2 = 0;
    double min_added_dist2 = 1e18;
    
    for (int i = 0; i <= new_hp.trips[trip2].drops.size(); i++) {
        Point prev = (i == 0)
            ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 ]
            : problem.villages[new_hp.trips[trip2].drops[i - 1].village_id].coords;
            
        Point next = (i == new_hp.trips[trip2].drops.size())
            ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 ]
            : problem.villages[new_hp.trips[trip2].drops[i].village_id].coords;
            
        double added_dist = distance(prev, drop1_point) + distance(drop1_point, next) - distance(prev, next);
        if (added_dist < min_added_dist2) {
            min_added_dist2 = added_dist;
            best_pos2 = i;
        }
    }
    
    // Insert the swapped drops
    new_hp.trips[trip1].drops.insert(new_hp.trips[trip1].drops.begin() + best_pos1, drop2);
    new_hp.trips[trip2].drops.insert(new_hp.trips[trip2].drops.begin() + best_pos2, drop1);
    
    // Update supplies with new drops
    new_hp.trips[trip1].dry_food_pickup += drop2.dry_food;
    new_hp.trips[trip1].perishable_food_pickup += drop2.perishable_food;
    new_hp.trips[trip1].other_supplies_pickup += drop2.other_supplies;
    
    new_hp.trips[trip2].dry_food_pickup += drop1.dry_food;
    new_hp.trips[trip2].perishable_food_pickup += drop1.perishable_food;
    new_hp.trips[trip2].other_supplies_pickup += drop1.other_supplies;
    
    if (!isFeasible(problem, new_hp)) {
        cout << "Swap resulted in infeasible plan" << endl;
        return hp;
    }
    
    return new_hp;
}

HelicopterPlan getRemoveNeighbor(const ProblemData& problem, const HelicopterPlan& hp){
    int helicopter_id = hp.helicopter_id;
    vector<Trip> trips = hp.trips;
    
    if (trips.empty()) return hp;
    
    random_device rd;
    mt19937 gen(rd());
    
    // Find trips with drops
    vector<int> valid_trips;
    for (int i = 0; i < trips.size(); i++) {
        if (!trips[i].drops.empty()) {
            valid_trips.push_back(i);
        }
    }
    
    if (valid_trips.empty()) {
        cout << "No trips with drops to remove" << endl;
        return hp;
    }
    
    // Select random trip
    uniform_int_distribution<> trip_dis(0, valid_trips.size() - 1);
    int trip_idx = valid_trips[trip_dis(gen)];
    
    // Select worst drop (highest distance from previous point)
    int worst_drop_idx = -1;
    double worst_score = -1e18;
    Point prev = problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1];
    
    for (int i = 0; i < trips[trip_idx].drops.size(); i++) {
        int village_id = trips[trip_idx].drops[i].village_id;
        Point curr = problem.villages[village_id].coords;
        double score = distance(curr, prev);
        
        if (score > worst_score) {
            worst_drop_idx = i;
            worst_score = score;
        }
        prev = curr;
    }
    
    if (worst_drop_idx == -1) {
        cout << "Could not find drop to remove" << endl;
        return hp;
    }
    
    HelicopterPlan new_hp = hp;
    Drop removed_drop = new_hp.trips[trip_idx].drops[worst_drop_idx];
    
    // Remove the drop
    new_hp.trips[trip_idx].drops.erase(new_hp.trips[trip_idx].drops.begin() + worst_drop_idx);
    
    // Update supplies
    new_hp.trips[trip_idx].dry_food_pickup -= removed_drop.dry_food;
    new_hp.trips[trip_idx].perishable_food_pickup -= removed_drop.perishable_food;
    new_hp.trips[trip_idx].other_supplies_pickup -= removed_drop.other_supplies;
    
    // Remove empty trips
    if (new_hp.trips[trip_idx].drops.empty()) {
        new_hp.trips.erase(new_hp.trips.begin() + trip_idx);
    }
    
    return new_hp;
}

HelicopterPlan getAddNeighbor(const ProblemData& problem, const HelicopterPlan& hp, vector<VillageAlt>& villages_alt) {
    int helicopter_id = hp.helicopter_id;
    
    // Find villages not currently assigned to this helicopter
    set<int> assigned_villages;
    for (const auto& trip : hp.trips) {
        for (const auto& drop : trip.drops) {
            assigned_villages.insert(drop.village_id);
        }
    }
    
    vector<int> unassigned_villages;
    for (int i = 0; i < problem.villages.size(); i++) {
        int village_id = problem.villages[i].id;
        if (assigned_villages.find(village_id) == assigned_villages.end()) {
            unassigned_villages.push_back(village_id);
        }
    }
    
    if (unassigned_villages.empty()) {
        cout << "No unassigned villages to add" << endl;
        return hp;
    }
    
    random_device rd;
    mt19937 gen(rd());
    
    // Select random unassigned village
    uniform_int_distribution<> village_dis(0, unassigned_villages.size() - 1);
    int selected_village_id = unassigned_villages[village_dis(gen)];
    
    // Create a drop for this village (you might need to adjust this based on your Drop structure)
    Drop new_drop;
    new_drop.village_id = selected_village_id;
    // Set supply amounts - you might want to use village requirements or default values
    const Village& village = problem.villages[selected_village_id - 1]; // assuming 1-indexed
    // find the alt village
    for (auto& alt_village : villages_alt) {
        if (alt_village.village.id == selected_village_id) {
            new_drop.dry_food = alt_village.food_required() / 2;
            new_drop.perishable_food = alt_village.food_required() / 2;
            new_drop.other_supplies = alt_village.other_required();
            break;
        }
    }

    HelicopterPlan new_hp = hp;
    
    // If no trips exist, create a new trip
    if (new_hp.trips.empty()) {
        Trip new_trip;
        new_trip.drops.push_back(new_drop);
        new_trip.dry_food_pickup = new_drop.dry_food;
        new_trip.perishable_food_pickup = new_drop.perishable_food;
        new_trip.other_supplies_pickup = new_drop.other_supplies;
        new_hp.trips.push_back(new_trip);
    } else {
        // Find best trip and position to add the drop
        int best_trip = -1;
        int best_pos = -1;
        double min_added_dist = 1e18;
        Point drop_point = problem.villages[selected_village_id - 1].coords;
        
        for (int t = 0; t < new_hp.trips.size(); t++) {
            // Try adding to existing trip
            for (int i = 0; i <= new_hp.trips[t].drops.size(); i++) {
                Point prev = (i == 0)
                    ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1]
                    : problem.villages[new_hp.trips[t].drops[i - 1].village_id - 1].coords;
                    
                Point next = (i == new_hp.trips[t].drops.size())
                    ? problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 ]
                    : problem.villages[new_hp.trips[t].drops[i].village_id - 1].coords;
                    
                double added_dist = distance(prev, drop_point) + distance(drop_point, next) - distance(prev, next);
                if (added_dist < min_added_dist) {
                    min_added_dist = added_dist;
                    best_trip = t;
                    best_pos = i;
                }
            }
        }
        
        // Add to best position
        if (best_trip != -1) {
            new_hp.trips[best_trip].drops.insert(new_hp.trips[best_trip].drops.begin() + best_pos, new_drop);
            new_hp.trips[best_trip].dry_food_pickup += new_drop.dry_food;
            new_hp.trips[best_trip].perishable_food_pickup += new_drop.perishable_food;
            new_hp.trips[best_trip].other_supplies_pickup += new_drop.other_supplies;
        }
    }
    
    if (!isFeasible(problem, new_hp)) {
        cout << "Adding drop resulted in infeasible plan" << endl;
        return hp;
    }
    
    return new_hp;
}

HelicopterPlan get2OptEfficientNeighbor(const ProblemData& problem, const HelicopterPlan& hp) {
    int helicopter_id = hp.helicopter_id;
    vector<Trip> trips = hp.trips;
    
    if (trips.empty()) return hp;
    
    random_device rd;
    mt19937 gen(rd());
    
    vector<int> valid_trips;
    for (int i = 0; i < trips.size(); i++) {
        if (trips[i].drops.size() >= 3) {
            valid_trips.push_back(i);
        }
    }
    
    if (valid_trips.empty()) return hp;
    
    uniform_int_distribution<> trip_dis(0, valid_trips.size() - 1);
    int selected_trip = valid_trips[trip_dis(gen)];
    
    vector<Drop>& drops = trips[selected_trip].drops;
    int n = drops.size();
    
    Point home = problem.cities[problem.helicopters[helicopter_id - 1].home_city_id - 1 ];
    

    double best_improvement = 0.0;
    int best_i = -1, best_j = -1;
    
    for (int i = 0; i < n - 2; i++) {
        for (int j = i + 2; j < n; j++) {
            
            
            Point pi = (i == 0) ? home : problem.villages[drops[i-1].village_id - 1].coords;
            Point pi1 = problem.villages[drops[i].village_id - 1].coords;
            Point pj = problem.villages[drops[j].village_id - 1].coords;
            Point pj1 = (j == n-1) ? home : problem.villages[drops[j+1].village_id - 1].coords;
            
            // Old edges: pi -> pi+1 and pj -> pj+1
            double old_dist = distance(pi, pi1) + distance(pj, pj1);
            
            // New edges: pi -> pj and pi+1 -> pj+1
            double new_dist = distance(pi, pj) + distance(pi1, pj1);
            
            double improvement = old_dist - new_dist;
            
            if (improvement > best_improvement) {
                best_improvement = improvement;
                best_i = i;
                best_j = j;
            }
        }
    }
    
    if (best_improvement > 0.001) { // Small epsilon to avoid floating point issues
        HelicopterPlan new_hp = hp;
        reverse(new_hp.trips[selected_trip].drops.begin() + best_i + 1, 
                new_hp.trips[selected_trip].drops.begin() + best_j + 1);
        
        cout << "Efficient 2-opt improvement: saved " << best_improvement << endl;
        return new_hp;
    }
    
    cout << "No efficient 2-opt improvement found" << endl;
    return hp;
}
/**
 * @brief The main function to implement your search/optimization algorithm.
 * * This is a placeholder implementation. It creates a simple, likely invalid,
 * plan to demonstrate how to build the Solution object. 
 * * TODO: REPLACE THIS ENTIRE FUNCTION WITH YOUR ALGORITHM.
 */

vector<HelicopterPlan> localSearch(const ProblemData& problem, Solution initial_plans, vector<VillageAlt>& villages_alt, 
    int max_iterations = 20) {
    vector<HelicopterPlan> current_plans = initial_plans;
    double current_score = evaluateSolution(problem, current_plans); 
    cout << "current_score = " << current_score << endl;
    double initial_temp = 100.0;
    double cooling_rate = 0.9;
    for (int iter = 0; iter < max_iterations; iter++) {

        for (int heli_idx = 0; heli_idx < current_plans.size(); heli_idx++) {

        // Generate a neighbor for that helicopter
        HelicopterPlan neighbor_plan1 = getCutPasteNeighbor(problem, current_plans[heli_idx]);
        HelicopterPlan neighbor_plan2 = getSwapNeighbor(problem, current_plans[heli_idx]);
        HelicopterPlan neighbor_plan3 = getRemoveNeighbor(problem, current_plans[heli_idx]);
        HelicopterPlan neighbor_plan4 = getAddNeighbor(problem, current_plans[heli_idx], villages_alt);
        HelicopterPlan neighbor_plan5 = get2OptEfficientNeighbor(problem, current_plans[heli_idx]);

        vector<HelicopterPlan> candidate_plans = current_plans;
        
        // see if any neighbor is strictly greater, if not then accept each neighbor with probability proportional to their exponential score
        double best_neighbor_score = evaluateHelicopterPlan(problem, current_plans[heli_idx], current_plans);
        HelicopterPlan best_neighbor;
        bool flag = 0;

        for (const HelicopterPlan& neighbor : {neighbor_plan1, neighbor_plan2, neighbor_plan3, neighbor_plan4, neighbor_plan5}) {
            double neighbor_score = evaluateHelicopterPlan(problem, neighbor, current_plans);
            if (neighbor_score > best_neighbor_score) {
                best_neighbor_score = neighbor_score;
                best_neighbor = neighbor;
                cout << "better neighbor found with score: " << neighbor_score << endl;
                flag = 1;
            }
        }
        if (flag) {
            current_plans[heli_idx] = best_neighbor;
        }
        // else {
        //     // select among the worse neighbors with a probability distribution
        //     // select among the worse neighbors with a probability distribution
        //     vector<HelicopterPlan> worse_neighbors {};
        //     for (const HelicopterPlan& neighbor : {neighbor_plan1, neighbor_plan2, neighbor_plan3, neighbor_plan4, neighbor_plan5}) {
        //         if (evaluatePlan(problem, neighbor) < best_neighbor_score) {
        //             worse_neighbors.push_back(neighbor);
        //             cout << "there are worse neighbors\n";
        //             cout << "Neighbor score: " << evaluatePlan(problem, neighbor) << endl;
        //         }
        //     }

        //     if (!worse_neighbors.empty()) {
        //         // If there's only one worse neighbor, select it directly
        //         if (worse_neighbors.size() == 1) {
        //             current_plans[heli_idx] = worse_neighbors[0];
        //             cout << "worse neighbor selected with score: " << evaluatePlan(problem, worse_neighbors[0]) << endl;
        //         } else {
        //             // Sample from multiple worse neighbors based on their scores, use temperatures next
        //             vector<double> scores;
        //             for (auto& nbr: worse_neighbors) {
        //                 double nbr_score = exp(evaluatePlan(problem, nbr) - evaluatePlan(problem, current_plans[heli_idx]));
        //                 scores.push_back(nbr_score);
        //             }
        //             // normalize scores
        //             double total_score = accumulate(scores.begin(), scores.end(), 0.0);
        //             for (auto& score : scores) {
        //                 score /= total_score;
        //             }
        //             // select the score index with probability equal to the score value
        //             discrete_distribution<> dist(scores.begin(), scores.end());
        //             int selected_index = dist(gen);
        //             current_plans[heli_idx] = worse_neighbors[selected_index];
        //             cout << "selected worse neighbor with score: " << evaluatePlan(problem, worse_neighbors[selected_index]) << endl;

        //         }
        //     }
        // }
        current_score = evaluateSolution(problem, current_plans);
        
        cout << "Iteration " << iter << ": score = " << current_score << endl;
        cout << "climbed: " << flag << endl;
    }
}

    return current_plans;
}
Solution solve(const ProblemData& problem) {

    auto create_priority_queue = [&](vector<VillageAlt>& villages_alt, Point curr, unordered_set<int>& visited_villages_ids) {
        priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, PQByHeuristic> pq;

        for(VillageAlt& v : villages_alt) {
            if(visited_villages_ids.count(v.village.id)) continue;

            double dist = distance(curr, v.village.coords);
            double heuristic;
            int food_required = v.food_required();
            int other_required = v.other_required();
            
            // double food_division_ratio = problem.packages[0].value / (problem.packages[0].value 
            // + problem.packages[1].value); // ratio of dry to total food in heuristic
            double food_division_ratio = 0.5; 
            heuristic = -dist;
            pq.push(make_pair(&v, heuristic));
        }

        return pq;
    };

    auto start = chrono::high_resolution_clock::now();

    cout << "Starting solver..." << endl;

    Solution solution;

    // sort(helicopters_alt.begin(), helicopters_alt.end(), [](const HelicopterAlt& a, const HelicopterAlt& b) {
    //     return a.helicopter.distance_capacity > b.helicopter.distance_capacity;
    // });

    double best_score = 0.0;

    unordered_map<int, HelicopterPlan> best_solution_map;

    for(Helicopter h : problem.helicopters){
            HelicopterPlan hp = {h.id, vector<Trip>()};
            best_solution_map[h.id] = hp;
        }

    unordered_map<string, vector<HelicopterAlt>> seen;

    int repeat_count = 0;
    vector<VillageAlt> villages_alt;

    while(time_check(start, problem.time_limit_minutes) && repeat_count < 100) {
        cout << "Time remaining: " << problem.time_limit_minutes - chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start).count() << " minutes" << endl;

        villages_alt.clear();
        for(Village v : problem.villages){
            VillageAlt va = {v, 0, 0, 0, false};
            villages_alt.push_back(va);
        }

        vector<HelicopterAlt> helicopters_alt;

        for(Helicopter h : problem.helicopters){
            HelicopterAlt ha = {h, 0.0, 0.0};
            helicopters_alt.push_back(ha);
        }

        //random shuffle
        random_device rd;
        mt19937 g(rd());
        shuffle(helicopters_alt.begin(), helicopters_alt.end(), g);

        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        // shuffle(helicopters_alt.begin(), helicopters_alt.end(), default_random_engine(seed));

        string key = vecToString(helicopters_alt);

        if (seen.find(key) == seen.end()) {
            // Store only new permutations
            seen[key] = helicopters_alt;
            repeat_count = 0;
        }
        else{
            repeat_count++;
            continue;
        }


        unordered_map<int, HelicopterPlan> solution_map;

        for(Helicopter h : problem.helicopters){
            HelicopterPlan hp = {h.id, vector<Trip>()};
            solution_map[h.id] = hp;
        }

        double score = 0.0;

        while(!satisfaction(villages_alt) && time_check(start, problem.time_limit_minutes)){
            int empty_drops_count = 0;
            int empty_trip_count = 0;

            double cumulative_trip_cost = 0.0;
            double cumulative_value_achieved = 0.0;
            
            for(HelicopterAlt& h : helicopters_alt){
                vector<Drop> drops;
                int dry_food_pickup = 0;
                int perishable_food_pickup = 0;
                int other_supplies_pickup = 0;

                double curr_weight = 0;

                Point curr = problem.cities[h.helicopter.home_city_id - 1];

                // vector<VillageAlt> visited_villages;
                unordered_set<int> visited_villages_ids;

                int village_check_count = 0;

                double trip_cost = 0.0;
                double value_achieved = 0.0;

                double prev_cumulative_distance = h.cumulative_distance_travelled;

                while(h.cumulative_distance_travelled < problem.d_max && h.trip_distance < h.helicopter.distance_capacity && village_check_count < villages_alt.size()){
                    priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, PQByHeuristic> pq = create_priority_queue(villages_alt, curr, visited_villages_ids);
                    if(pq.empty()) break;
                    while(!pq.empty()){
                        VillageAlt& v = *pq.top().first;
                        pq.pop();

                        village_check_count++;

                        if(v.is_fulfilled) continue;
                        if(distance(curr, v.village.coords) + distance(v.village.coords, problem.cities[h.helicopter.home_city_id - 1 ]) + h.trip_distance > h.helicopter.distance_capacity) continue;
                        if(distance(curr, v.village.coords) + distance(v.village.coords, problem.cities[h.helicopter.home_city_id - 1 ]) + h.cumulative_distance_travelled > problem.d_max) continue;

                        //check for requirements of village
                        int food_required = v.village.population * 9 - (v.dry_food_recieved + v.perishable_food_recieved);
                        int other_required = v.village.population - v.other_supplies_recieved;
                        double remaining_capacity = h.helicopter.weight_capacity - curr_weight;

                        //greedy wrt type of food
                        //first perishable food
                        int perishable_food_add = min(food_required, (int) (remaining_capacity / problem.packages[1].weight));
                        perishable_food_pickup += perishable_food_add;
                        curr_weight += perishable_food_add * problem.packages[1].weight;
                        food_required -= perishable_food_add;
                        remaining_capacity = h.helicopter.weight_capacity - curr_weight;

                        //second dry food
                        int dry_food_add = min(food_required, (int) (remaining_capacity/ problem.packages[0].weight));
                        dry_food_pickup += dry_food_add;
                        curr_weight += dry_food_add * problem.packages[0].weight;
                        food_required -= dry_food_add;
                        remaining_capacity = h.helicopter.weight_capacity - curr_weight;

                        //third other supplies
                        int other_supplies_add = min(other_required, (int) (remaining_capacity / problem.packages[2].weight));
                        other_supplies_pickup += other_supplies_add;
                        curr_weight += other_supplies_add * problem.packages[2].weight;
                        other_required -= other_supplies_add;
                        remaining_capacity = h.helicopter.weight_capacity - curr_weight;

                        v.dry_food_recieved += dry_food_add;
                        v.perishable_food_recieved += perishable_food_add;
                        v.other_supplies_recieved += other_supplies_add;
                        // bool f = false; // changed
                        if(dry_food_add + perishable_food_add + other_supplies_add > 0){
                            h.cumulative_distance_travelled += distance(curr, v.village.coords);
                            h.trip_distance += distance(curr, v.village.coords);
            
                            curr = v.village.coords;
                            visited_villages_ids.insert(v.village.id);

                            Drop drop = {v.village.id, dry_food_add, perishable_food_add, other_supplies_add};
                            drops.push_back(drop);
                            // f = true; // changed
                        }

                        bool food_ok = (v.dry_food_recieved + v.perishable_food_recieved) >= 9 * v.village.population;
                        bool other_ok = v.other_supplies_recieved >= v.village.population;
                        if(food_ok && other_ok) v.is_fulfilled = true;
                        // if (f) break; // changed
                    }
                }

                if(dry_food_pickup + perishable_food_pickup + other_supplies_pickup == 0){
                    empty_trip_count++;
                    continue;
                }

                h.cumulative_distance_travelled += distance(curr, problem.cities[h.helicopter.home_city_id - 1]);
                h.trip_distance = 0;

                trip_cost += h.helicopter.fixed_cost + h.helicopter.alpha * (h.cumulative_distance_travelled - prev_cumulative_distance);
                value_achieved += dry_food_pickup * problem.packages[0].value + perishable_food_pickup * problem.packages[1].value + other_supplies_pickup * problem.packages[2].value;

                if(drops.size() == 0){
                    empty_drops_count ++;
                }
                else{
                    Trip trip = {dry_food_pickup, perishable_food_pickup, other_supplies_pickup, drops};
                    solution_map[h.helicopter.id].trips.push_back(trip);

                    cout << "Heli" << h.helicopter.id << " " << dry_food_pickup << " " << perishable_food_pickup << " " << other_supplies_pickup << " " << h.cumulative_distance_travelled << endl;
                }

                cumulative_trip_cost += trip_cost;
                cumulative_value_achieved += value_achieved;
            }

            score += cumulative_value_achieved - cumulative_trip_cost;

            //printing output
            cout << "Village details" << endl;
            for(VillageAlt v : villages_alt){
                cout << v.village.id << " " << v.dry_food_recieved << " " << v.perishable_food_recieved << " " << v.other_supplies_recieved << endl;
            }
            //printing end

            if(empty_drops_count == helicopters_alt.size()){
                //no more drops can be made
                break;
            }

            if(empty_trip_count == helicopters_alt.size()){
                //no more trips can be made
                break;
            }
        }

        if(score > best_score){
            best_score = score;
            best_solution_map = solution_map;

            cout << "New best score: " << best_score << endl;
        }
    }

    for(Helicopter h : problem.helicopters){
        solution.push_back(best_solution_map[h.id]);
    }

    

    cout << "Solver finished." << endl;
    // return solution;
    cout << best_score << "*****" << endl;

    Solution final_solution {solution};
    Solution local_search_solution = localSearch(problem, solution, villages_alt, 70);
    double score_after_local_search = evaluateSolution(problem, local_search_solution);
    best_score = evaluateSolution(problem, final_solution);
    cout << best_score << "******" << endl;
    if (score_after_local_search > best_score) {
        cout << "New best score after local search: " << score_after_local_search << endl;
        return local_search_solution;
    } else {
        cout << "No improvement found after local search." << endl;
        return final_solution;
    }

}

