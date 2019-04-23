#include <iostream>
#include "tbb/task_scheduler_init.h"
#include "tbb/task_group.h"
#include "tbb/mutex.h"
#include "tbb/tick_count.h"
using namespace tbb;

#define Eps 0.01



double val_one_dimensional_func(double x);
double val_two_dimensonal_func(double x1, double x2);
double val_two_dimensonal_func_2(double x1, double x2);

void one_dimensional_integral(double a, double b, double h, double *res);
void one_dim_integral_parallel(double a, double b, double h, double *res);

void two_dimensional_integral_fix_point(double a1, double a2, double b2, double h, double *res);
void two_dimensional_integral_multi_step(double a1, double b1, double a2, double b2, double h, double *res);

void CheckResultsTest(double sequentTime, double parallelTime);

class TaskOneDimIntegral : public task {
 private:
     double left_;
     double right_;
     double h_;
     double *res_;
     double local_res_;
 public:
     TaskOneDimIntegral(double left, double right, double h, double *res) : left_(left), right_(right), h_(h), res_(res) {} 
     task* execute() {
         if (right_ - left_ < 20.0) {
             local_res_ = 0.0;
             one_dimensional_integral(left_, right_, h_, &local_res_);
             *res_ += local_res_;
             return NULL;
         }
         else {
             double local_res1 = 0.0, local_res2 = 0.0;
             task_list list; 
             TaskOneDimIntegral& a = *new(task::allocate_child()) TaskOneDimIntegral(left_, left_ +  (right_ - left_) / 2, h_, &local_res1);
             TaskOneDimIntegral& b = *new(task::allocate_child()) TaskOneDimIntegral(left_ + (right_ - left_) / 2, right_, h_, &local_res2);
             list.push_back(a);
             list.push_back(b);
             task::set_ref_count(3);
             spawn_and_wait_for_all(list);
             *res_ += local_res1;
             *res_ += local_res2;
         }
         return NULL;
     }

     double get_res() {
         return *res_;
     }
};

class TaskTwoDimIntegral : public task {
 private:
    double x_left_;
    double x_right_;
    double y_left_;
    double y_right_;
    double h_;
    double *res_;

 public:
     TaskTwoDimIntegral(double x_left, double x_right, double y_left, double y_right,
         double h, double *res) : x_left_(x_left), x_right_(x_right), y_left_(y_left), y_right_(y_right),h_(h), res_(res) {}

     task* execute() {
         if (x_right_ - x_left_ < 1.0) {
             double local_res = 0.0;
             two_dimensional_integral_multi_step(x_left_, x_right_, y_left_,
                 y_right_,h_, &local_res);
             *res_ += local_res;
             return NULL;
         }
         else {
             double local_res1 = 0.0, local_res2 = 0.0;
             task_list list;
             TaskTwoDimIntegral& a = *new(task::allocate_child()) TaskTwoDimIntegral(x_left_, x_left_ + (x_right_ - 
                 x_left_) / 2, y_left_, y_right_, h_, &local_res1);
             TaskTwoDimIntegral& b = *new(task::allocate_child()) TaskTwoDimIntegral(x_left_ + (x_right_ 
                 - x_left_) / 2, x_right_, y_left_, y_right_, h_, &local_res2);
             list.push_back(a);
             list.push_back(b);
             task::set_ref_count(3);
             spawn_and_wait_for_all(list);
             *res_ += local_res1;
             *res_ += local_res2;
         }
         return NULL;
     }
};

class TaskTwoDimIntegralXY : public task {
private:
    double x_left_;
    double x_right_;
    double y_left_;
    double y_right_;
    double h_;
    double *res_;
    double local_res_;

public:
    TaskTwoDimIntegralXY(double x_left, double x_right, double y_left, double y_right,
        double h, double *res) : x_left_(x_left), x_right_(x_right), y_left_(y_left), 
        y_right_(y_right), h_(h), res_(res) {}

    task* execute() {
        if (x_right_ - x_left_ <= 1.0 && y_right_ - y_left_ <= 1.0) {
            local_res_ = 0.0;
            two_dimensional_integral_multi_step(x_left_, x_right_, y_left_,
                y_right_, h_, &local_res_);
            *res_ += local_res_;
            return NULL;
        }
        else if (x_right_ - x_left_ <= 1.0 && y_right_ - y_left_ > 1.0) {
            double local_res1 = 0.0, local_res2 = 0.0;
            task_list list;
            TaskTwoDimIntegralXY& a = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_, x_right_,
                y_left_, y_left_ + (y_right_ - y_left_) / 2, h_, &local_res1);
            TaskTwoDimIntegralXY& b = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_, x_right_,
                y_left_ + (y_right_ - y_left_) / 2, y_right_, h_, &local_res2);
            list.push_back(a);
            list.push_back(b);
            task::set_ref_count(3);
            spawn_and_wait_for_all(list);
            *res_ += local_res1;
            *res_ += local_res2;
        }
        else if (x_right_ - x_left_ > 1.0 && y_right_ - y_left_ <= 1.0) {
            double local_res1 = 0.0, local_res2 = 0.0;
            task_list list;
            TaskTwoDimIntegralXY& a = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_, x_left_ +
                (x_right_ - x_left_) / 2, y_left_, y_right_, h_, &local_res1);
            TaskTwoDimIntegralXY& b = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_ +
                (x_right_ - x_left_) / 2, x_right_, y_left_, y_right_, h_, &local_res2);
            list.push_back(a);
            list.push_back(b);
            task::set_ref_count(3);
            spawn_and_wait_for_all(list);
            *res_ += local_res1;
            *res_ += local_res2;
        }
        else {
            double local_res1 = 0.0, local_res2 = 0.0, local_res3 = 0.0, local_res4 = 0.0;
            task_list list;
            TaskTwoDimIntegralXY& a = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_, x_left_ + (x_right_ -
                x_left_) / 2, y_left_, y_left_ + (y_right_ - y_left_) / 2, h_, &local_res1);
            TaskTwoDimIntegralXY& b = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_ + (x_right_
                - x_left_) / 2, x_right_, y_left_, y_left_ + (y_right_ - y_left_) / 2, h_, &local_res2);
            TaskTwoDimIntegralXY& c = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_, x_left_ + (x_right_ -
                x_left_) / 2, y_left_ + (y_right_ - y_left_) / 2, y_right_, h_, &local_res3);
            TaskTwoDimIntegralXY& d = *new(task::allocate_child()) TaskTwoDimIntegralXY(x_left_ + (x_right_
                - x_left_) / 2, x_right_, y_left_ + (y_right_ - y_left_) / 2, y_right_, h_, &local_res4);

            list.push_back(a);
            list.push_back(b);
            list.push_back(c);
            list.push_back(d);
            task::set_ref_count(5);
            spawn_and_wait_for_all(list);
            *res_ += local_res1;
            *res_ += local_res2;
            *res_ += local_res3;
            *res_ += local_res4;
        }
        return NULL;
    }
};



double val_one_dimensional_func(double x) {
    return sin(x);
}
double val_two_dimensonal_func(double x1, double x2) {
    return sin(x1) * cos(x2);
}
double val_two_dimensonal_func_2(double x1, double x2) {
    return 0.5*sin(x1*x2)*x1*x2*log((x1/10) + 1);
}

void one_dimensional_integral(double a, double b, double h, double *res) {
    double sum = 0.0;
    for (double i = a; i < b; i += h) {
        sum += val_one_dimensional_func(i) * h;
    }
    *res = sum;
}
void one_dim_integral_parallel(double a, double b, double h, double *res) {
    double total_res = 0.0;
    int i;
    for (i = 0; i < int((b - a) / h); i++) {
        total_res += val_one_dimensional_func(a + i * h) * h;
    }

    *res = total_res;
}


void two_dimensional_integral_multi_step(double a1, double b1, double a2, double b2, double h, double *res) {
    double sum = 0.0;
    double tmpsum = 0.0;
    int i;
    for (i = 0; i < int((b1 - a1) / h); i++)
    {
        two_dimensional_integral_fix_point(a1 + h * i, a2, b2, h, &tmpsum);
        sum += tmpsum * h;
    }
    *res = sum;
}

void two_dimensional_integral_fix_point(double a1, double a2, double b2, double h, double *res) {
    double sum = 0.0;

    for (double j = a2; j < b2; j += h)
    {
        sum += val_two_dimensonal_func_2(a1, j) * h;
    }

    *res = sum;
}

void CheckResultsTest(double sequentTime, double parallelTime) {
    std::cout << "The analysis of the running time of algorithms " << std::endl;
    if (parallelTime < sequentTime)
        std::cout << "Parallel algorithm is faster" << std::endl;
    else
        std::cout << "Sequential algorithm is faster" << std::endl;

    std::cout << "Speedup " << sequentTime / parallelTime << std::endl;
}

int main(int argc, char **argv)
{

    double a1 = 0.0, b1 = 10.0, a2 = 0.0, b2 = 10.0;
    double h = 0.001;
    double seqRes = 0.0, parRes = 0.0;
    int num_threads = 8;
    double result = 0.0;
    double seqTime = 0.0, parTime = 0.0;

    tbb::task_scheduler_init init(num_threads);
    tick_count t0;
    tick_count t1;

    for (int i = 1; i < argc; ++i) {
        if ((!strcmp(argv[i], "-a1")) && (i + 1 < argc)) {
            a1 = atof(argv[i + 1]);
        }
        if ((!strcmp(argv[i], "-b1")) && (i + 1 < argc)) {
            b1 = atof(argv[i + 1]);
        }
        if ((!strcmp(argv[i], "-h")) && (i + 1 < argc)) {
            h = atof(argv[i + 1]);
        }
        if ((!strcmp(argv[i], "-a2")) && (i + 1 < argc)) {
            a2 = atof(argv[i + 1]);
        }
        if ((!strcmp(argv[i], "-b2")) && (i + 1 < argc)) {
            b2 = atof(argv[i + 1]);
        }
    }

    if (a2 == INFINITY || b2 == INFINITY) {
        std::cout << "*****calculating the value of the one-dimensional integral*****" << std::endl;
        t0 = tick_count::now();
        one_dimensional_integral(a1, b1, h, &seqRes);
        t1 = tick_count::now();
        std::cout << "Work took sequential version " << (t1 - t0).seconds() << " seconds" << std::endl;
        t0 = tick_count::now();
        TaskOneDimIntegral& task = *new(task::allocate_root()) TaskOneDimIntegral(a1, b1, h, &parRes);
        task::spawn_root_and_wait(task);
        t1 = tick_count::now();
        std::cout << "Work took parallel version " << (t1 - t0).seconds() << " seconds" << std::endl;
        /*for (int i = 0; i < num_threads; i++) {
            size++;
            left = a1 + i * step;
            right = left + step;
            TaskOneDimIntegral& task = *new(task::allocate_root()) TaskOneDimIntegral(left, right, h, &local_res);
            tasks.push_back(task);
        }
        task::spawn_root_and_wait(tasks);
        result = local_res;*/
        /*for (int i = 0; i < num_threads; i++) {
            left = a1 + i * step;
            right = left + step;
            local_res = 0;
            TaskOneDimIntegral& task = *new(task::allocate_root()) TaskOneDimIntegral(left, right, h, local_res);
            task::spawn_root_and_wait(task);
            //tg.run(TaskOneDimIntegral(left, right, h, local_res));
            std::cout << "Local Result " << local_res << std::endl;
            result += local_res;
        }
        tg.wait();*/
        std::cout << "Result one dimersion integral (Parallel version) " << parRes << std::endl;

        //CheckResultsTest(sequentTimeWork, parallTimeWork);
        //Test(seqRes, parRes);
    }
    else {
        std::cout << "*****calculating the value of the two-dimensional integral*****" << std::endl;

        t0 = tick_count::now();
        two_dimensional_integral_multi_step(a1, b1, a2, b2, h, &seqRes);
        t1 = tick_count::now();
        seqTime = (t1 - t0).seconds();
        std::cout << "Work took sequential version " << seqTime << " seconds" << std::endl;
        std::cout << "Result one dimersion integral (sequential version) " << seqRes << std::endl;

        t0 = tick_count::now();
        TaskTwoDimIntegral& task = *new(task::allocate_root()) TaskTwoDimIntegral(a1, b1, a2, b2, h, &parRes);
        task::spawn_root_and_wait(task);
        t1 = tick_count::now();
        
        parTime = (t1 - t0).seconds();
        std::cout << "Work took parallel version " << parTime << " seconds" << std::endl;
        std::cout << "Result two dimersion integral (Parallel version) " << parRes << std::endl << std::endl;
        parRes = 0.0;

        CheckResultsTest(seqTime, parTime);

        std::cout << "*********** Dividing the interval by X and Y ************" << std::endl;
        

        t0 = tick_count::now();
        TaskTwoDimIntegralXY& taskNew = *new(task::allocate_root()) TaskTwoDimIntegralXY(a1, b1, a2, b2, h, &parRes);
        task::spawn_root_and_wait(taskNew);
        t1 = tick_count::now();
        std::cout << "Work took parallel version " << (t1 - t0).seconds() << " seconds" << std::endl;
        std::cout << "Result two dimersion integral (Parallel version) " << parRes << std::endl;

        //Test(seqRes, parRes);
    }


}