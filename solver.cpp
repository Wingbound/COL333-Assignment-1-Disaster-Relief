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
};

struct HelicopterAlt{
    Helicopter helicopter;
    double cumulative_distance_travelled;
    double trip_distance;
};

struct PQByDist {
    bool operator()(const pair<VillageAlt*, double>& a,
                    const pair<VillageAlt*, double>& b) const {
        return a.second > b.second; // min-heap by distance
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
    return elapsed.count() < time_limit_minutes * 60 * 1000 - 10*1000;
}

priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, PQByDist> unvisited_villages_queue(vector<VillageAlt>& villages_alt, Point curr, unordered_set<int>& visited_villages_ids){
    priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, PQByDist> pq;

    for(VillageAlt& v : villages_alt){
        if(visited_villages_ids.count(v.village.id)) continue;

        double dist = distance(curr, v.village.coords);
        pq.push(make_pair(&v, dist));
    }

    return pq;
}

string vecToString(vector<HelicopterAlt>& helicopters_permutations) {
    string s;
    for (HelicopterAlt h : helicopters_permutations) {
        s += to_string(h.helicopter.id);
    }
    return s;
}


/**
 * @brief The main function to implement your search/optimization algorithm.
 * * This is a placeholder implementation. It creates a simple, likely invalid,
 * plan to demonstrate how to build the Solution object. 
 * * TODO: REPLACE THIS ENTIRE FUNCTION WITH YOUR ALGORITHM.
 */
Solution solve(const ProblemData& problem) {

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

    while(time_check(start, problem.time_limit_minutes) && repeat_count < 100){
        cout << "Time remaining: " << problem.time_limit_minutes - chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start).count() << " minutes" << endl;

        vector<VillageAlt> villages_alt;

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
                    priority_queue<pair<VillageAlt*, double>, vector<pair<VillageAlt*, double>>, PQByDist> pq = unvisited_villages_queue(villages_alt, curr, visited_villages_ids);
                    if(pq.empty()) break;
                    while(!pq.empty()){
                        VillageAlt& v = *pq.top().first;
                        pq.pop();

                        village_check_count++;

                        if(v.is_fulfilled) continue;
                        if(distance(curr, v.village.coords) + distance(v.village.coords, problem.cities[h.helicopter.home_city_id - 1]) + h.trip_distance > h.helicopter.distance_capacity) continue;
                        if(distance(curr, v.village.coords) + distance(v.village.coords, problem.cities[h.helicopter.home_city_id - 1]) + h.cumulative_distance_travelled > problem.d_max) continue;

                        //check for requirements of village
                        int food_required = v.village.population*9 - (v.dry_food_recieved + v.perishable_food_recieved);
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

                        if(dry_food_add + perishable_food_add + other_supplies_add > 0){
                            h.cumulative_distance_travelled += distance(curr, v.village.coords);
                            h.trip_distance += distance(curr, v.village.coords);
            
                            curr = v.village.coords;
                            visited_villages_ids.insert(v.village.id);

                            Drop drop = {v.village.id, dry_food_add, perishable_food_add, other_supplies_add};
                            drops.push_back(drop);
                        }

                        bool food_ok = (v.dry_food_recieved + v.perishable_food_recieved) >= 9 * v.village.population;
                        bool other_ok = v.other_supplies_recieved >= v.village.population;
                        if(food_ok && other_ok) v.is_fulfilled = true;

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
    cout << "Final score: " << best_score << endl;

    cout << "Solver finished." << endl;
    return solution;
}