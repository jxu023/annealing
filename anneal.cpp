#include <math.h>
#include <set>
#include <queue>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <climits>
#include <list>

using namespace std;

const bool debug = true;
struct SimulatedAnnealing {
    // generic params read in
    double initprob;
    int sizefactor;
    double cutoff;
    double tempfactor;
    int freeze_lim;
    double minpercent;
    // vars
    const int CHROM_EST = 10; // increasing increases runtime
    int m,n; // edges, nodes
    int neighborhood_size;
    double c_star;
    vector<int> solution_star;
    int colors_used;
    double c;
    vector<int> solution;
    double c_prime;
    vector<int> solution_prime;
    int FIXED_K;
    bool stochastic_tunneling;
    double gamma;
    double lowestcost;
    int iteration;
    // data structures
    vector<list<int>> adj;
    map<int,int> color_size;
    map<int,int> color_bad_edges;
    map<int,int> color_size_prime;
    map<int,int> color_bad_edges_prime;
    int approach;
    // problem specific subroutines
    void read_instance(string filename) {
        ifstream ifh(filename);
        string line;
        string tmp;
        char c;
        while (getline(ifh,line)) {
            if (line[0] == 'p' || line[0] == 'e'){
                stringstream ss(line);
                if (line[0] == 'p') {
                    ss >> c; ss >> tmp;
                    ss >> n; ss >> m;
                    adj.resize(n,list<int>());
                    solution_star.resize(n,0);
                    solution.resize(n,0);
                    solution_prime.resize(n,0);
                }
                else {
                    ss >> c;
                    int v1, v2;
                    ss >> v1; ss >> v2;
                    adj[v1-1].push_back(v2-1);
                    adj[v2-1].push_back(v1-1);
                }
            }
        }
        ifh.close();
        c_star = INT_MAX; // todo : set to upper bound
        c_prime = INT_MAX;
        colors_used = -1;
        neighborhood_size = n * CHROM_EST;
        if (approach == 3)
            neighborhood_size = n * FIXED_K;
    }
    bool islegal;
    double cost(vector<int> & sol) {
        color_size_prime.clear();
        double cost1 = 0.0;
        for (int cl : sol)
            ++color_size_prime[cl];
        if (approach != 3) {
            for (auto iter : color_size_prime) {
                auto value = iter.second;
                cost1 -= (double)(value)*(value);
            }
        }
        color_bad_edges_prime.clear();
        islegal = true;
        if (approach == 1 || approach == 3) {
            for (int v1 = 0; v1 < n; ++v1) {
                for (int v2 : adj[v1]) {
                    if (sol[v1] == sol[v2]) {
                        ++color_bad_edges_prime[sol[v1]];
                        if (approach == 3) {
                            cost1 += 1.0;
                        }
                        islegal = false;
                    }
                }
            }
            if (approach == 3) {
                if (cost1 > 0.0)
                    cost1 = cost1 / 2.0;
            }
            else {
                for (auto iter : color_size_prime) {
                    // the coeff of 2 is already inside the count of bad edges ... since ea edge counted twice
                    cost1 += (double)iter.second*color_bad_edges_prime[iter.first];
                }
            }
        }
        if (stochastic_tunneling) {
            if (cost1 != 0.0)
                cost1 *= (1.0 - exp(-1.0*fabs(cost1-lowestcost)/gamma));
        }
        return cost1;
    }
    void change_soln() {
        c = c_prime;
        solution = solution_prime;
        color_size = color_size_prime;
        color_bad_edges = color_bad_edges_prime;
        if (islegal && c_prime < c_star) {
            ofstream ofh("anneal.trace",std::ios_base::app);
            ofh << iteration << ' ' << c_prime << '\n';
            ofh.close();
            changed = true;
            solution_star = solution_prime;
            c_star = c_prime;
            colors_used = color_size.size();
        }
    }
    void initial_solution() {
        if (approach == 3) {
            for (int i = 0; i < n; ++i) {
                solution[i] = rand() % FIXED_K;
            }
        } else {
            for (int i = 0; i < n; ++i) {
                solution[i] = i;
            }
        }
        /*
        for (int i = 0; i < n; ++i) {
            solution[i] = rand() % CHROM_EST;
        }
        c = cost(solution);
        // do color shifting to make nonempty colors continguous
        map<int,int> m;
        int i = 0;
        for (auto iter : color_size_prime) {
            m[iter.first] = i++;
        }
        for (int & num : solution) {
            num = m[num];
        }
        */
        c = cost(solution);
        color_size = color_size_prime;
        color_bad_edges = color_bad_edges_prime;

        solution_prime = solution;
        c_prime = c;
        lowestcost = c;

        change_soln();

        ofstream ofh("anneal.trace");
        ofh << c << '\n';
        ofh.close(); // clear the trace
    }
    // penalty function swap
    void change1() {
        int k = color_size.size();
        int color = rand() % (k+1);
        int vert = rand() % n;
        solution_prime = solution;
        int old_color = solution_prime[vert];
        while (old_color == color || (color_size[old_color] == 1 && color == k+1)) {
            color = rand() % (k+1);
        }
        solution_prime[vert] = color;
        // shift colors
        if (color_size[old_color] == 1) {
            for (auto & num : solution_prime)
                if (num > old_color)
                    --num;
        }
        c_prime = cost(solution_prime);
    }
    // kempe chain swap
    void change2() {
        map<int,bool> full_chain; // take product of two vertices to index this
        int k = color_size.size();
        int num_full = k*(k-1)/2; // number of possible pairs
        set<int> kempe_chain; // set is used to manage duplicates
        int color1, color2;
        solution_prime = solution;
        while(1) {
            kempe_chain.clear();
            int vert = rand() % n; 
            color1 = solution[vert];
            while ((color2=rand()%k) == color1 || full_chain[(color1+1)*(color2+1)]);
            unsigned union_size = color_size[color1] + color_size[color2];
            kempe_chain.insert(vert);
            queue<int> q;
            q.push(vert);
            while(!q.empty()) {
                vert = q.front();
                q.pop();
                for (int v2 : adj[vert]) {
                    if ((solution[v2] == color2 || solution[v2] == color1) && kempe_chain.find(v2) == kempe_chain.end()){
                        kempe_chain.insert(v2);
                        q.push(v2);
                    }
                }
            }
            if (kempe_chain.size() < union_size)
                break;
            full_chain[(color1+1)*(color2+1)] = true;
            if (--num_full <= 0) {
                cout << "there are no more possible kempe chain swaps" << endl;
                change_soln();
                final_soln();
                exit(0);
            }
        }
        list<int> assignC;
        int assignedC = 0;
        for (int v1 : kempe_chain) {
            if (solution[v1] == color2) {
                assignC.push_back(v1);
                ++assignedC;
            }
            else {
                solution_prime[v1] = color2;
                --assignedC;
            }
        }
        for (int v1 : assignC)
            solution_prime[v1] = color1;
        // check for empty color classes, do color index shift
        if (color_size[color1] <= -1*assignedC || color_size[color2] <= assignedC) {
            int oldcolor = color2;
            if (color_size[color1] <= -1*assignedC) // note that it should never be <, only == or >
                oldcolor = color1;
            for (int & num : solution_prime)
                if (num > oldcolor)
                    --num;
        }
        /*
        cout << "solprime: ";
        for (int num : solution_prime) {
            cout << num << ' ';
        }
        cout << endl;*/
        c_prime = cost(solution_prime);
    }
    void change3() {
        solution_prime = solution;
        set<int> bad_vert;
        for (int v1 = 0; v1 < n; ++v1) {
            for (int v2 : adj[v1]) {
                if (solution[v1] == solution[v2]) {
                    bad_vert.insert(v1);
                    bad_vert.insert(v2);
                }
            }
        }
        int B = bad_vert.size();
        int ind = rand() % B;
        auto it = bad_vert.begin();
        advance(it,ind);
        int vert = *it;
        int oldcolor = solution[vert];
        int color = rand() % (FIXED_K-1);
        if (color >= oldcolor)
            ++color;
        solution_prime[vert] = color;
        c_prime = cost(solution_prime); // update should be simple for this.. oh well
    }
    void next_change() {
        if (approach == 1) {
            change1();
        }
        else if (approach == 2) {
            change2();
        }
        else if (approach == 3) {
            change3();
        }
        lowestcost = min(lowestcost,c_prime);
    }
    bool changed;
    void final_soln() {
        ofstream ofh("anneal.run");
        ofh << c_star << '\n';
        ofh << colors_used << '\n';
        for (int i = 0; i < n; ++i) {
            ofh << i << ' ' << solution_star[i] << '\n';
        }
        ofh.close();
        if (debug) {
            cout << "final_soln()\n";
            cout << "iteration: " << iteration << endl;
            cout << "solution : \n";
            for (int num : solution) {
                cout << num << ' ';
            }
            cout << "solution star: \n";
            for (int num : solution_star) {
                cout << num << ' ';
            }
            cout << endl;
            cout << "best cost is: " << c_star << endl;
            if (approach == 3) {
                cout << "FIXED_K is " << FIXED_K << endl;
            }
            cout << "colors used: " << colors_used << endl;

            if (c_star == INT_MAX)
                cout << "\n\nNO LEGAL SOL FOUND\n\n";
        }
    }
    double find_T() {
        double l = 0.001; // could make this lower
        double r = 1000;
        cout << "\nsearching for T right now\n\n";
        iteration=0;
        while (l < r) {
            double mid = (l+r)/2;
            double T = mid;
            initial_solution();
            if (approach == 3 && islegal) { // **
                cout << "FIXED K initial solution is legal" << endl;
                return 0.0;
            } // ***
            int changes, trials;
            changes = trials = 0;
            int count = 0;
            while (trials < sizefactor*neighborhood_size && changes < cutoff * neighborhood_size && count < 1000) {
                ++trials;
                ++iteration;
                next_change();
                if (approach == 3 && islegal) {
                    if (debug)
                        cout << "found fixed k sol in trial run\n";
                    return 0;
                }
                double delta = c_prime - c;
                if (delta <= 0) {
                    ++changes;
                    change_soln();
                }
                else {
                    double r = (double)rand() / RAND_MAX;
                    if (r <= exp(-1.0*delta/T)) {
                        ++changes;
                        change_soln();
                    }
                }
                ++count;
            }
            if (debug) {
                cout << "trials: " << trials << " changes: " << changes << endl;
                cout << "prob is " << (double)changes/trials << endl;
                cout << "trial find T, current T is " << T << endl;
            }
            double val = (double)changes/trials;
            if (val > 0.5 || val == 0.0) {
                l = 10.0;
                break;
            }
            const double diff = 0.001;
            if (val < initprob)
                l = mid + diff;
            else
                r = mid - diff;
        }
        cout << "\n\npicked T value of " << l << endl;
        return l;
    }
    // generic simulated annealing
    void run_sa(string filename) {
        cout << fixed;
        cout << setprecision(10);
        iteration = 0;
        read_instance(filename);
        double T = find_T();
        if (T == 0.0) {
            change_soln();
            final_soln();
            return;
        }
        initial_solution();
        if (approach == 3 && islegal) { // **
            cout << "FIXED K initial solution is legal" << endl;
            change_soln();
            final_soln();
            return;
        } // ***
        int freeze_count = 0;
        int changes,trials;
        vector<int> prev(n,-1);
        int same_changes = 0;
        int same_cost = 0;
        int prevchanges = -1;
        double prevcost = INT_MAX;
        while(freeze_count < freeze_lim) {
            changes = trials = 0;
            changed = false;
            while (trials < sizefactor*neighborhood_size && changes < cutoff * neighborhood_size) {
                ++trials;
                ++iteration;
                next_change();
                if (approach == 3 && (islegal || iteration >= 200000)) {
                    change_soln();
                    final_soln();
                    return;
                }
                double delta = c_prime - c;
                if (delta <= 0) {
                    ++changes;
                    change_soln();
                }
                else {
                    double r = (double)rand() / RAND_MAX;
                    if (r <= exp(-1.0*delta/T)) {
                        ++changes;
                        change_soln();
                    }
                }
            }
            if (changes == prevchanges)
                ++same_changes;
            else {
                prevchanges = changes;
                same_changes = 0;
            }
            if (c == prevcost)
                ++same_cost;
            else {
                prevcost = c;
                same_cost = 0;
            }
            if (same_changes > 50 || same_cost > 50) { // prevents inf looping for kempe chain
                // maybe dynamic stoch tunneling would solve this
                break;
            }
            T = tempfactor*T;
            if ((double)changes / trials < minpercent)
                ++freeze_count;
            if (changed)
                freeze_count = 0;
            if (debug) {
                cout << "\n\ntrials " << trials << "\n";
                cout << "iteration " << iteration << "\n";
                cout << "current cost is " << c << endl;
                cout << "current num colors is  " << color_size.size() << endl;
                cout << "best cost is " << c_star << endl;
                cout << "lowest colors is " << colors_used << endl;
                cout << "changes is " << changes << endl;
                cout << "freeze_count is " << freeze_count << endl;
                cout << "same_changes is " << same_changes << endl;
                cout << "same_cost is " << same_cost << endl;
            }
            if (approach == 3) { // prevents inf loops in some cases
                int i;
                for (i = 0; i < n; ++i) {
                    if (solution[i] != prev[i])
                        break;
                }
                if (i == n)
                    break;
                prev = solution;
            }
        }
        final_soln();
    }
    void set_params() {
        ifstream ifh("params.in");
        string name;
        while (ifh.good()) {
            ifh >> name;
            if (name == "INITPROB") {
                ifh >> initprob;
            }
            else if (name == "FREEZE_LIM") {
                ifh >> freeze_lim;
            }
            else if (name == "SIZEFACTOR") {
                ifh >> sizefactor;
            }
            else if (name == "CUTOFF") {
                ifh >> cutoff;
            }
            else if (name == "TEMPFACTOR") {
                ifh >> tempfactor;
            }
            else if (name == "MINPERCENT") {
                ifh >> minpercent;
            }
            else if (name == "GAMMA") {
                ifh >> gamma;
            }
        }
        ifh.close();
    }
    /*
    ASSUMPTION:
    undirected graph, input file and thus adj list has ea edge appear just once (e.g. e 1 2, no e 2 1)

    TODOS:
    todos for organization
        put the data structures/vars associated with a single solution into a struct of its own
            this includes cost value, list of bad edges per color, # of nodes per color ...
            implement structure history functionality, allowing for "undos"
            this avoids expensive copy operations that are clearly inefficient for simple perturbations
        multiple files
            create a virtual class SA, and 3 classes that inherit from SA to impl ea approach
                virtual class SA contains stoch tunneling option
                this improves efficiency too
        create a testing harness to run multiple tests at once

    todos for better correctness
        ********************
        Use a map<pair<int,int>,bool> for FULL CHAIN detection .. otherwise 3,4 and 2,6 store to the same thing
        or use a map<int,set<int>>, v1 is key, if v2 is in the returned set then return true
        *******************
        impl. DSATUR or rand seq coloring for initial solution
            this ensures that there exists a value of T causing changes/trials to be close to INITPROB
            if the initial sol. is really bad, most of the neighboring sol is good
            this would cause changes/trials to be unable to go below some value say 0.5 no matter what T is
            for a good solution with mostly bad neighbors, we can adjust T to allow more/less bad solutions as changes
            then we can set T with much greater precision
        avoid "relooping" while picking colors, just reduce the random range and "map up"

    todos for efficiency
        add problem specific update_error() func ... avoid iterating over every edge for ea perturb/trial
            auto-differentiation for update error func?
                fancy .. uneccessary unless generalizing this stuff much more
        use one vector instead of two vectors for solution and solution_prime
            just need to have an undo operation
            this avoids the much more expensive copy operation in next_change()
        remove the branches in the functions (next_change, init_sol)
            by adding function pointers
        tabu search
            this can be really low cost and easy if just store a hash of the previous 100 solutions
            decide not to move to a solution if its hash is the same as a prev one
        read instance is called more times than necessary for fixed_k
    */
    SimulatedAnnealing(string inputfile, int seed, int appr, bool st,int fixed) {
        approach = appr;
        stochastic_tunneling = st;
        srand(seed);
        set_params();
        FIXED_K = fixed;
        run_sa(inputfile);
    }
};

int main(int argc, char ** argv)
{
    string filename;
    int seed;
    int approach;
    bool st=false;
    if (argc == 4 || argc == 5) {
        filename = argv[1];
        seed = atoi(argv[2]);
        string appr = argv[3];
        if (appr == "-penalty") {
            approach = 1;
        }
        else if (appr == "-kempe") {
            approach = 2;
        }
        else
            approach = 3;
        if (argc == 5)
            st = true;
        if (approach == 3) {
            int l = 2;
            int r = 100;
            while (l < r) {
                int mid = (l+r)/2;
                cout << "testing k of " << mid << endl;
                SimulatedAnnealing sa(filename,seed,approach,st,mid);
                if (sa.islegal) {
                    r = mid;
                }
                else {
                    l = mid+1;
                }
            }
            cout << "min k is " << r << endl;
        }
        else
            SimulatedAnnealing sa(filename,seed,approach,st,0);
    }
    else {
        cout << "invalid args! run as:\n./anneal <filename> <seed> {-penalty,-kempe,-fixed_k} [-st]";
    }
    return 0;
}
