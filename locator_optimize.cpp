#include <stdio.h>

#include "ibs_matrix.h"
#include "ga/ga.h"
#include "ga/std_stream.h"
#include "ga/GAGenome.h"
#include "ga/gaversion.h"
#include "mex_function.h"
#include "mex-it/mex-it.h"
#include <stack>
#include <algorithm>
#include <python2.7/Python.h>

using namespace std;

const double north[] = {1, 0};

double dot_product(const point_t &veca, const point_t &vecb) {
    return veca.x * vecb.x + veca.y * vecb.y;
}

double norm(const point_t &vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y);
}

point_t cross(const homo_point_t &from, const homo_point_t &to) {
    double x = from.y * to.z - from.z * to.y;
    double y = from.z * to.x - from.x * to.z;
    double z = from.x * to.y - from.y * to.x;

    return point_t(x / z, y / z);
}

double pdist(const point_t &veca, const point_t &vecb) {
    return sqrt((veca.x - vecb.x) * (veca.x - vecb.x) +
                (veca.y - vecb.y) * (veca.y - vecb.y));
}

point_t p0(0, 0);

bool compares(const point_t &p1, const point_t &p2) {
    double t;
    t = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
    if (t > 0 || (t == 0 && pdist(p1, p0) < pdist(p2, p0)))
        return true;
    else
        return false;
}

double area(point_t p1, point_t p2, point_t p3) {
    double s;
    s = (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)) / 2.0;
    return s;
}

// remove the kth point
double get_area_from_point_set(vector<point_t> input, int k) {
    stack<point> s;

    if (k >= 0) {
        input.erase(input.begin() + k);
    }

    int n = input.size();

    double x = input[0].x;
    double y = input[0].y;
    int p = 0;
    for (int i = 1; i < n; i++)      //寻找最左上的点
    {
        if ((input[i].y == y && input[i].x < x) || input[i].y < y) {
            x = input[i].x;
            y = input[i].y;
            p = i;
        }
    }
    //x,y为最做点，即p0确定
    p0.x = x;
    p0.y = y;

    swap(input[0], input[p]);
    sort(input.begin() + 1, input.end(), compares);
    //将0、1、2 三个点压入堆栈 S
    s.push(p0);
    s.push(input[1]);
    s.push(input[2]);

    double p_temp;
    point_t temp(0, 0);
    for (int i = 3; i < n; i++) {       //凸包
        while (!s.empty()) {
            temp = s.top();
            s.pop();
            p_temp = (s.top().x) * temp.y + input[i].x * (s.top().y) + temp.x * (input[i].y) - input[i].x * (temp.y) -
                     temp.x * (s.top().y) - s.top().x * (input[i].y);
            if (p_temp > 0) { //i位于spot[t-1]spot[t]左侧
                s.push(temp); //将temp重新压入栈中
                break;
            }
        }
        s.push(input[i]); //将i 压入栈中
    }

    double area = 0;
    temp = s.top();
    s.pop();
    area += (p0.x * temp.y - p0.y * temp.x) / 2.0;
    while (!s.empty()) {
        area += (temp.x * s.top().y - temp.y * s.top().x) / 2.0;
        temp = s.top();
        s.pop();
    }
    if (n <= 2) {
        area = 0;
    }
    if (area < 0) {
        area = -area;
    }

    return area;
}

/*****************************
 ** begin genetic algorithm **
 *****************************/

float objective(GAGenome &);

vector_d genome2vector(GAGenome &);

vector_d calculate_location_result(const vector_d &delta_theta);

int landmark_count;
vector<point_t> landmark_location;
vector_d compass;

vector_d
ga_main(const vector_d &match_result_x, const vector_d &match_result_y, const vector_d &compass_reading) {

    printf("------GA VERSION: %s ------\n", GA_VERSION);

    landmark_location.clear();
    compass.clear();

    assert(match_result_x.size() == match_result_y.size());
    assert(match_result_x.size() == compass_reading.size());
    landmark_count = match_result_x.size();

    printf("Landmark number is: %d\n", landmark_count);
    printf("Location data: \n");
    for (int i = 0; i < landmark_count; i++) {
        point_t p(match_result_x[i], match_result_y[i]);
        landmark_location.push_back(p);
        compass.push_back(compass_reading[i]);
        printf("%f\t%f\t%f\n", p.x, p.y, compass[i]);
    }

// See if we've been given a seed to use (for testing purposes).  When you
// specify a random seed, the evolution will be exactly the same each time
// you use that seed number.

    unsigned int seed = time(NULL);

// Declare variables for the GA parameters and set them to some default values.

    int popsize = 30;
    int ngen = 200;

    float pmut = 0.08;
    float pcross = 0.9;

// Create a phenotype for two variables.  The number of bits you can use to
// represent any number is limited by the type of computer you are using.

    GABin2DecPhenotype map;
    for (int i = 0; i < landmark_count; ++i) {
        map.add(16, GA_MIN_VALUE, GA_MAX_VALUE);
    }

// Create the template genome using the phenotype map we just made.

    GABin2DecGenome genome(map, objective);

// Now create the GA using the genome and run it.  We'll use sigma truncation
// scaling so that we can handle negative objective scores.

    GASimpleGA ga(genome);
    GASigmaTruncationScaling scaling;
    ga.minimize();
    ga.populationSize(popsize);
    ga.nGenerations(ngen);
    ga.pMutation(pmut);
    ga.pCrossover(pcross);
    ga.scaling(scaling);
    ga.scoreFilename("locator.dat");
    ga.scoreFrequency(1);
    ga.flushFrequency(10);
    ga.evolve(seed);

// Dump the results of the GA to the screen.
    while (!ga.done()) {
        ga.step();
    }

    genome = ga.statistics().bestIndividual();
    cout << "the ga found an optimum at the point (";
    for (int j = 0; j < genome.nPhenotypes() - 1; ++j) {
        cout << genome.phenotype(j) << ", ";
    }
    cout << genome.phenotype(genome.nPhenotypes() - 1) << ")\n\n";

    vector_d delta = genome2vector(genome);
    vector_d location_result = calculate_location_result(delta);

    printf("Location Result is:(%f, %f), cost is: %f, %f, %f, %f\n", location_result[0], location_result[1],
           location_result[2], location_result[3], location_result[4], location_result[5]);
    return location_result;
}

vector_d genome2vector(GAGenome &ga) {
    GABin2DecGenome &genome = (GABin2DecGenome &) ga;

    vector_d delta;
    for (int i = 0; i < genome.nPhenotypes(); ++i) {
        delta.push_back(genome.phenotype(i));
    }
    return delta;
}

int counter = 0;

void lof_from_python(const vector<point_t> &input, vector<point_t> output) {
    PyObject *pModule = PyImport_Import(PyString_FromString("lof"));

    PyObject *input_list = PyList_New(input.size());

    for (int i = 0; i < input.size(); i++) {
        PyObject *tuple = PyTuple_New(2);
        PyTuple_SetItem(tuple, 0, PyFloat_FromDouble(input[i].x));
        PyTuple_SetItem(tuple, 1, PyFloat_FromDouble(input[i].y));
        PyList_SetItem(input_list, i, tuple);
    }

    PyObject *pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, PyInt_FromLong(input.size()));
    PyTuple_SetItem(pArgs, 1, input_list);

//    printf("Input size is: %d, output size is: %d\n", PyList_Size(input_list), PyInt_AsLong(PyTuple_GetItem(pArgs, 0)));

    PyObject *pFunc = PyObject_GetAttrString(pModule, "outliers");

//    printf("Callable: %d\n", PyCallable_Check(pFunc));

    PyObject *pValue = PyObject_CallObject(pFunc, pArgs);

    printf("Input size is: %d, output size is: %d\n", PyList_Size(input_list), PyList_Size(pValue));
}


vector_d calculate_location_result(const vector_d &delta_theta) {
    vector<point_t> cross_result;
    vector<double> weight;

    for (int i = 0; i < landmark_count - 1; ++i) {
        for (int j = i + 1; j < landmark_count; ++j) {
            point_t landmark_a = landmark_location[i];
            point_t landmark_b = landmark_location[j];

            double theta_a = compass[i] + delta_theta[i];
            double theta_b = compass[j] + delta_theta[j];

            iBS::Matrix rotation_a(2), rotation_b(2);
            rotation_a.data[0][0] = cos(DEG2RAD(theta_a));
            rotation_a.data[0][1] = sin(DEG2RAD(theta_a));
            rotation_a.data[1][0] = -sin(DEG2RAD(theta_a));
            rotation_a.data[1][1] = cos(DEG2RAD(theta_a));

            rotation_b.data[0][0] = cos(DEG2RAD(theta_b));
            rotation_b.data[0][1] = sin(DEG2RAD(theta_b));
            rotation_b.data[1][0] = -sin(DEG2RAD(theta_b));
            rotation_b.data[1][1] = cos(DEG2RAD(theta_b));

            iBS::Matrix north_vec(2);
            north_vec.data[0][0] = north[0];
            north_vec.data[1][0] = north[1];

            iBS::Matrix direction_vec_a = rotation_a * north_vec;
            iBS::Matrix direction_vec_b = rotation_b * north_vec;

            point_t dir_a(direction_vec_a.data[0][0], direction_vec_a.data[1][0]);
            point_t dir_b(direction_vec_b.data[0][0], direction_vec_b.data[1][0]);

            homo_point_t line_a(dir_a.y, -dir_a.x, dir_a.x * landmark_a.y - dir_a.y * landmark_a.x);
            homo_point_t line_b(dir_b.y, -dir_b.x, dir_b.x * landmark_b.y - dir_b.y * landmark_b.x);

            point_t user_pos = cross(line_a, line_b);

//            printf("landmark %d and landmark %d : (%f, %f)\n", i, j, user_pos.x, user_pos.y);

            // weight = (sind(theta(i) - theta(j)) .^ 4) / ...
            //             (pdist(store_pos([i, j], :))) .^ 2;

            cross_result.push_back(user_pos);
            weight.push_back(1.0 / (pdist(landmark_a, landmark_b) * pdist(landmark_a, landmark_b)));
//            weight.push_back(1.0 / (pdist(landmark_a, landmark_b)));
        }
    }

//    vector<point_t> filter_list;
//    lof_from_python(cross_result, filter_list);
//    printf("Count: %d\n", counter++);

    if (cross_result.size() >= 6) {
//        // calculate cross distance
        vector<vector<double> > distance(cross_result.size(), vector<double>(cross_result.size()));
        for (int i = 0; i < cross_result.size() - 1; i++) {
            for (int j = i + 1; j < cross_result.size(); j++) {
                distance[i][j] = distance[j][i] = pdist(cross_result[i], cross_result[j]);
            }
        }

        vector<double> sum_dist(cross_result.size());
        for (int i = 0; i < cross_result.size(); i++) {
            for (int j = 0; j < cross_result.size(); j++) {
                sum_dist[i] += distance[i][j];
            }
        }

        vector<int> index(cross_result.size());
        for (int i = 0; i < cross_result.size(); i++) { index[i] = i; }
        sort(index.begin(), index.end(), [&](int x, int y) { return sum_dist[x] > sum_dist[y]; });

        for (int i = 0; i < landmark_count - 1; i++) {
            if(sum_dist[index[i]] > sum_dist[index[i+1]] * 1.2) {
                weight[index[i]] = 0;
            }
        }
        // calculate distance
//        double total_area = get_area_from_point_set(cross_result, -1);
//        vector<double> area_difference(cross_result.size());
//
//        for (int i = 0; i < cross_result.size(); i++) {
//            area_difference[i] = total_area - get_area_from_point_set(cross_result, i);
//        }
//        vector<int> index(cross_result.size());
//        for (int i = 0; i < cross_result.size(); i++) { index[i] = i; }
//        sort(index.begin(), index.end(), [&](int x, int y) { return area_difference[x] > area_difference[y]; });
//
//        weight[index[0]] = 0;
    }

    // normalize weight
    double weight_sum = 0;
    for (int i = 0; i < weight.size(); i++) {
        weight_sum += weight[i];
    }

    // calucuate average location result
    double average_x = 0;
    double average_y = 0;
    for (int i = 0; i < cross_result.size(); i++) {
        average_x += cross_result[i].x * (weight[i] / weight_sum);
        average_y += cross_result[i].y * (weight[i] / weight_sum);
    }

    point_t loc_result(average_x, average_y);


    double loc_cost = 0;
    for (int i = 0; i < cross_result.size(); ++i) {
        loc_cost += pdist(loc_result, cross_result[i])
                    * pdist(loc_result, cross_result[i]) * (weight[i] / weight_sum);
//        loc_cost += pdist(loc_result, cross_result[i]) * (weight[i] / weight_sum);
    }
//    loc_cost /= cross_result.size();

//    cout << "the ga search at point: (";
//    for (int i = 0; i < delta_theta.size() - 1; i++) {
//        cout << delta_theta[i] << ", ";
//    }
//    cout << delta_theta[delta_theta.size() - 1] << ")" << endl;

//    cout << "loc_cost: " << loc_cost << endl;

    double delta_cost = 0;
    for (int i = 0; i < delta_theta.size(); i++) {
        delta_cost += delta_theta[i] * delta_theta[i];
//        delta_cost += abs(delta_theta[i]);
    }
//    cout << "delta_cost: " << delta_cost << endl;

    double landmark_cost = 0;
    for (int i = 0; i < landmark_location.size(); i++) {
        landmark_cost += pdist(loc_result, landmark_location[i]) * pdist(loc_result, landmark_location[i]);
//        landmark_cost += pdist(loc_result, landmark_location[i]);
    }
    landmark_cost /= landmark_location.size();
//    cout << "landmark_cost: " << landmark_cost << endl;

    double total_cost = LAMBDA_A * loc_cost + LAMBDA_B * delta_cost + LAMBDA_C * landmark_cost;
//    cout << "total_cost: " << total_cost << endl;

    vector_d result;
    result.push_back(loc_result.x);
    result.push_back(loc_result.y);
    result.push_back(total_cost);
    result.push_back(loc_cost);
    result.push_back(delta_cost);
    result.push_back(landmark_cost);
    return result;
}

float
objective(GAGenome &c) {
    GABin2DecGenome &genome = (GABin2DecGenome &) c;

    vector_d delta = genome2vector(genome);

    vector_d result = calculate_location_result(delta);

    // result vector: loc_x, loc_y, cost;
    return result[2];
}


// We just need the mex_function header file
// This could be put into mex-it.h but this way allows naming of the mex function in matlab to this file name (without the .cxx)
void mex_function(const vector_d &match_result_x, const vector_d &match_result_y, const vector_d &compass_reading,
                  vector_d &location_result) {
//    Py_Initialize();
//
//    if (!Py_IsInitialized()) {
//        printf("ERROR!");
//        return;
//    }
//
//    PyRun_SimpleString("import sys");
//    PyRun_SimpleString("sys.path.append(\"/data/UserData/Mingkuan/SourceCode/0214-WALLE/pylof/pylof\")");
//    PyRun_SimpleString("print(sys.path)");

    location_result = ga_main(match_result_x, match_result_y, compass_reading);

//    Py_Finalize();
}

