//
//
//  Task_2.45_Comp._Prak.
//
//  Created by Denis Ashurov on 12.10.18.
//  Copyright © 2018 Denis Ashurov. All rights reserved.
//


#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define alpha 10
#define mu 1.0/3.5
#define la 4.0
#define facmax 2.0
#define facmin 1.e-10
#define fac 0.9
#define pi 3.1415926

double f_1(double x, double y, double t); //Both equations
double f_2(double x, double y, double t);
double min(double a, double b);
double max(double a, double b);
double runge_1(double x, double y, double t, double h); // Both functions for evaluation of new value of x(t), y(t)
double runge_2(double x, double y, double t, double h);
double solution(double x, double y, double t, double eps); // Equation solver for different initial conditions, returns x(t)
double solution_1(double x, double y, double t, double L, double eps);
double solution_2(double x, double y, double t, double L, double eps);// Clone, returning y(t)
double cyclepoint(void); // Cyclepoint searcher, returns cyclepoint
double final_solution(double x, double y, double t, double eps); //Getting data for plotting periodic solution
void runge_numbers(double x, double t ); // Function evaluating Runge numbers
double e_lambdamax(double x, double y, double h);

int main(int argc, const char * argv[]) {
    double k, t, T;
    t = 0;
    k = 0;
    k = cyclepoint();
    printf("Cyclepoint = %.10lf\n", k);
    t = final_solution(k, 0, t, 1.e-14) - t;
    printf("Period = %.9lf\n", t);
    T = t; //Period
    runge_numbers(k, 50*pi);
    return 0;
}

double f_1(double x, double y, double t) {
    return x + y - alpha*pow(x,3);
    t = t;
}

double f_2(double x, double y, double t) {
    return -x + y - pow(y,3);
    t = t;
}

double min(double a, double b) {
    if(a - b < 0) return a;
    else return b;
}

double max(double a, double b) {
    if(a - b > 0) return a;
    else return b;
}

double runge_1(double x, double y, double t, double h) {
    double k_11, k_12, k_21, k_22, k_31, k_32, k_41, k_42, k_51, k_52, k_61, k_62, k_71, k_72;
    double c_2, c_3, c_4, c_5, c_6, c_7;
    double a_21, a_31, a_32, a_41, a_42, a_43, a_51, a_52, a_53, a_54, a_61, a_62, a_63, a_64, a_65, a_71, a_72, a_73, a_74, a_75, a_76;
    double b_1, b_2, b_3, b_4, b_5, b_6, b_7; // объявление переменных матрицы Бутчера
    
    c_2 = mu;
    c_3 = 2.0/3.0;
    c_4 = 1.0/3.0;
    c_5 = 1.0/2.0;
    c_6 = 1.0/2.0;
    c_7 = 1;
    a_21 = mu;
    a_31 = 2.0/3.0 - 2.0/(9.0*mu);
    a_32 = 2.0/(9.0*mu);
    a_41 = 5.0/12.0 - 1.0/(9.0*mu);
    a_42 = 1.0/(9.0*mu);
    a_43 = -1.0/12.0;
    a_51 = 17.0/16.0 - 3.0/(8.0*mu);
    a_52 = 3.0/(8.0*mu);
    a_53 = - 3.0/16.0;
    a_54 = - 3.0/8.0;
    a_61 = 17.0/16.0 - 3.0/(8.0*mu) + 1.0/(4.0*la);
    a_62 = 3.0/(8.0*mu);
    a_63 = - 3.0/16.0 - 3.0/(4.0*la);
    a_64 = - 3.0/8.0 - 3.0/(2*la);
    a_65 = 2.0/la;
    a_71 = -27.0/44.0 + 3.0/(11.0*mu);
    a_72 = -3.0/(11.0*mu);
    a_73 = 63.0/44.0;
    a_74 = 18.0/11.0;
    a_75 = (4.0*la-16.0)/11.0;
    a_76 = - (4.0*la)/11.0;
    b_1 = 11.0/120.0;
    b_2 = 0;
    b_3 = 27.0/40.0;
    b_4 = 27.0/40.0;
    b_5 = (la - 8.0)/15.0;
    b_6 = -la/15.0;
    b_7 = 11.0/120.0;
    k_11 = f_1(x, y, t);
    k_12 = f_2(x, y, t);
    k_21 = f_1(x + a_21*h*k_11, y + a_21*h*k_12, t + c_2*h);
    k_22 = f_2(x + a_21*h*k_11, y + a_21*h*k_12, t + c_2*h);
    k_31 = f_1(x + a_31*h*k_11 + a_32*h*k_21, y + a_31*h*k_12 + a_32*h*k_22, t + c_3*h);
    k_32 = f_2(x + a_31*h*k_11 + a_32*h*k_21, y + a_31*h*k_12 + a_32*h*k_22, t + c_3*h);
    k_41 = f_1(x + a_41*h*k_11 + a_42*h*k_21 + a_43*h*k_31, y + a_41*h*k_12 + a_42*h*k_22 + a_43*h*k_32, t + c_4*h);
    k_42 = f_2(x + a_41*h*k_11 + a_42*h*k_21 + a_43*h*k_31, y + a_41*h*k_12 + a_42*h*k_22 + a_43*h*k_32, t + c_4*h);
    k_51 = f_1(x + a_51*h*k_11 + a_52*h*k_21 + a_53*h*k_31 + a_54*h*k_41, y + a_51*h*k_12 + a_52*h*k_22 + a_53*h*k_32 + a_54*h*k_42, t + c_5*h);
    k_52 = f_2(x + a_51*h*k_11 + a_52*h*k_21 + a_53*h*k_31 + a_54*h*k_41, y + a_51*h*k_12 + a_52*h*k_22 + a_53*h*k_32 + a_54*h*k_42, t + c_5*h);
    k_61 = f_1(x + a_61*h*k_11 + a_62*h*k_21 + a_63*h*k_31 + a_64*h*k_41 + a_65*h*k_51, y + a_61*h*k_12 + a_62*h*k_22 + a_63*h*k_32 + a_64*h*k_42 + a_65*h*k_52, t + c_6*h);
    k_62 = f_2(x + a_61*h*k_11 + a_62*h*k_21 + a_63*h*k_31 + a_64*h*k_41 + a_65*h*k_51, y + a_61*h*k_12 + a_62*h*k_22 + a_63*h*k_32 + a_64*h*k_42 + a_65*h*k_52, t + c_6*h);
    k_71 = f_1(x + a_71*h*k_11 + a_72*h*k_21 + a_73*h*k_31 + a_74*h*k_41 + a_75*h*k_51 + a_76*h*k_61, y + a_71*h*k_12 + a_72*h*k_22 + a_73*h*k_32 + a_74*h*k_42 + a_75*h*k_52 + a_76*h*k_62, t + c_7*h);
    k_72 = f_2(x + a_71*h*k_11 + a_72*h*k_21 + a_73*h*k_31 + a_74*h*k_41 + a_75*h*k_51 + a_76*h*k_61, y + a_71*h*k_12 + a_72*h*k_22 + a_73*h*k_32 + a_74*h*k_42 + a_75*h*k_52 + a_76*h*k_62, t + c_7*h);
    return x + h*(b_1*k_11 + b_2*k_21 + b_3*k_31 + b_4*k_41 + b_5*k_51 + b_6*k_61 + b_7*k_71);
}

double runge_2(double x, double y, double t, double h) {
    double k_11, k_12, k_21, k_22, k_31, k_32, k_41, k_42, k_51, k_52, k_61, k_62, k_71, k_72;
    double c_2, c_3, c_4, c_5, c_6, c_7;
    double a_21, a_31, a_32, a_41, a_42, a_43, a_51, a_52, a_53, a_54, a_61, a_62, a_63, a_64, a_65, a_71, a_72, a_73, a_74, a_75, a_76;
    double b_1, b_2, b_3, b_4, b_5, b_6, b_7; // объявление переменных матрицы Бутчера
    
    c_2 = mu;
    c_3 = 2.0/3.0;
    c_4 = 1.0/3.0;
    c_5 = 1.0/2.0;
    c_6 = 1.0/2.0;
    c_7 = 1;
    a_21 = mu;
    a_31 = 2.0/3.0 - 2.0/(9.0*mu);
    a_32 = 2.0/(9.0*mu);
    a_41 = 5.0/12.0 - 1.0/(9.0*mu);
    a_42 = 1.0/(9.0*mu);
    a_43 = -1.0/12.0;
    a_51 = 17.0/16.0 - 3.0/(8.0*mu);
    a_52 = 3.0/(8.0*mu);
    a_53 = - 3.0/16.0;
    a_54 = - 3.0/8.0;
    a_61 = 17.0/16.0 - 3.0/(8.0*mu) + 1.0/(4.0*la);
    a_62 = 3.0/(8.0*mu);
    a_63 = - 3.0/16.0 - 3.0/(4.0*la);
    a_64 = - 3.0/8.0 - 3.0/(2*la);
    a_65 = 2.0/la;
    a_71 = -27.0/44.0 + 3.0/(11.0*mu);
    a_72 = -3.0/(11.0*mu);
    a_73 = 63.0/44.0;
    a_74 = 18.0/11.0;
    a_75 = (4.0*la-16.0)/11.0;
    a_76 = - (4.0*la)/11.0;
    b_1 = 11.0/120.0;
    b_2 = 0;
    b_3 = 27.0/40.0;
    b_4 = 27.0/40.0;
    b_5 = (la - 8.0)/15.0;
    b_6 = -la/15.0;
    b_7 = 11.0/120.0;
    k_11 = f_1(x, y, t);
    k_12 = f_2(x, y, t);
    k_21 = f_1(x + a_21*h*k_11, y + a_21*h*k_12, t + c_2*h);
    k_22 = f_2(x + a_21*h*k_11, y + a_21*h*k_12, t + c_2*h);
    k_31 = f_1(x + a_31*h*k_11 + a_32*h*k_21, y + a_31*h*k_12 + a_32*h*k_22, t + c_3*h);
    k_32 = f_2(x + a_31*h*k_11 + a_32*h*k_21, y + a_31*h*k_12 + a_32*h*k_22, t + c_3*h);
    k_41 = f_1(x + a_41*h*k_11 + a_42*h*k_21 + a_43*h*k_31, y + a_41*h*k_12 + a_42*h*k_22 + a_43*h*k_32, t + c_4*h);
    k_42 = f_2(x + a_41*h*k_11 + a_42*h*k_21 + a_43*h*k_31, y + a_41*h*k_12 + a_42*h*k_22 + a_43*h*k_32, t + c_4*h);
    k_51 = f_1(x + a_51*h*k_11 + a_52*h*k_21 + a_53*h*k_31 + a_54*h*k_41, y + a_51*h*k_12 + a_52*h*k_22 + a_53*h*k_32 + a_54*h*k_42, t + c_5*h);
    k_52 = f_2(x + a_51*h*k_11 + a_52*h*k_21 + a_53*h*k_31 + a_54*h*k_41, y + a_51*h*k_12 + a_52*h*k_22 + a_53*h*k_32 + a_54*h*k_42, t + c_5*h);
    k_61 = f_1(x + a_61*h*k_11 + a_62*h*k_21 + a_63*h*k_31 + a_64*h*k_41 + a_65*h*k_51, y + a_61*h*k_12 + a_62*h*k_22 + a_63*h*k_32 + a_64*h*k_42 + a_65*h*k_52, t + c_6*h);
    k_62 = f_2(x + a_61*h*k_11 + a_62*h*k_21 + a_63*h*k_31 + a_64*h*k_41 + a_65*h*k_51, y + a_61*h*k_12 + a_62*h*k_22 + a_63*h*k_32 + a_64*h*k_42 + a_65*h*k_52, t + c_6*h);
    k_71 = f_1(x + a_71*h*k_11 + a_72*h*k_21 + a_73*h*k_31 + a_74*h*k_41 + a_75*h*k_51 + a_76*h*k_61, y + a_71*h*k_12 + a_72*h*k_22 + a_73*h*k_32 + a_74*h*k_42 + a_75*h*k_52 + a_76*h*k_62, t + c_7*h);
    k_72 = f_2(x + a_71*h*k_11 + a_72*h*k_21 + a_73*h*k_31 + a_74*h*k_41 + a_75*h*k_51 + a_76*h*k_61, y + a_71*h*k_12 + a_72*h*k_22 + a_73*h*k_32 + a_74*h*k_42 + a_75*h*k_52 + a_76*h*k_62, t + c_7*h);
    return y + h*(b_1*k_12 + b_2*k_22 + b_3*k_32 + b_4*k_42 + b_5*k_52 + b_6*k_62 + b_7*k_72);
}

double solution(double x, double y, double t, double eps) {
    double h, y_11, y_12, y_21, y_22, w_1, w_2, d_1, d_2, y_prev;
    bool flag = true; // flag
    double err;
    
    h = 0.3;  // first step
    y_11 = 0;
    y_12 = 0;
    y_21 = 0;
    y_22 = 0;
    
    y_prev = -1;
    
    while(flag) {
        y_11 = runge_1(x, y, t, h);
        y_12 = runge_2(x, y, t, h);
        y_21 = runge_1(y_11, y_12, t + h, h);
        y_22 = runge_2(y_11, y_12, t + h, h);
        w_1 = runge_1(x, y, t, 2*h);
        w_2 = runge_2(x, y, t, 2*h);
        d_1 = fabs(y_21 + (y_21 - w_1)/63);
        d_2 = fabs(y_22 + (y_22 - w_2)/63);
        err = max(fabs(y_21 - w_1)/d_1, fabs(y_22 - w_2)/d_2)/63.0;
        
        if( err < eps) {
            if(y_prev*y_22 < 0 && y_21 > 0) {
                h = h/2;
            }
            else {
                y = y_22;
                x = y_21;
                t = t + 2*h;
                h = h*min(facmax, max(facmin, fac*pow(eps/err, 1.0/7.0 )));
                y_prev = y_22;
            }
        }
        else {
            h = h*min(1.0, max(facmin, fac*pow(eps/err, 1.0/7.0)) );
        }
        if( fabs(y_22) < eps && y_21 > 0) flag = false;
    }
    
    return x;
}

double cyclepoint() {
    double x_start_1, x_start_2,t, y, x_end, cp;
    cp = 0;
    bool flag_2 = true;
    x_start_1 = 0.001;
    t = 0;
    y = 0;
    x_start_2 = solution(x_start_1, y, t, 1.e-8);
    x_end = solution(x_start_2, y, t, 1.e-8);
    while(flag_2) {
        if( fabs(x_end - x_start_2) < fabs(x_start_2 - x_start_1)) {
            x_start_1 = x_end;
            x_start_2 = solution(x_start_1, y, t, 1.e-8);
            x_end = solution(x_start_2, y, t, 1.e-8);
        }
        else {
            x_start_1 = x_start_1/2;
            x_start_2 = solution(x_start_1, y, t, 1.e-8);
            x_end = solution(x_start_2, y, t, 1.e-8);
        }
        if(fabs(x_end - x_start_2)< 1.e-12) flag_2 = false;
    }
    return x_end;
}

double final_solution(double x, double y, double t, double eps) {
    FILE *Out_x, *Out_y, *Out_t; // output files
    double h, y_11, y_12, y_21, y_22, w_1, w_2, d_1, d_2, y_prev;
    bool flag = true; // flag
    double err;
    
    h = 0.1;  // first step
    y_11 = 0;
    y_12 = 0;
    y_21 = 0;
    y_22 = 0;
    
    if(!(Out_x = fopen("output_x.txt", "w"))) {
        printf("Can't open output file for wrighting\n");
        return(-1);
    }
    
    if(!(Out_y = fopen("output_y.txt", "w"))) {
        printf("Can't open output file for wrighting\n");
        return(-1);
    }
    
    if(!(Out_t = fopen("output_t.txt", "w"))) {
        printf("Can't open output file for wrighting\n");
        return(-1);
    }
    
    fprintf(Out_x, "%5.11lf, ", x); // printing initial conditions
    fprintf(Out_y, "%5.11lf, ", y);
    fprintf(Out_t, "%5.11lf, ", t);
    y_prev = -1;
    
    while(flag) {
        y_11 = runge_1(x, y, t, h);
        y_12 = runge_2(x, y, t, h);
        y_21 = runge_1(y_11, y_12, t + h, h);
        y_22 = runge_2(y_11, y_12, t + h, h);
        w_1 = runge_1(x, y, t, 2*h);
        w_2 = runge_2(x, y, t, 2*h);
        d_1 = fabs(y_21 + (y_21 - w_1)/63);
        d_2 = fabs(y_22 + (y_22 - w_2)/63);
        err = max(fabs(y_21 - w_1)/d_1, fabs(y_22 - w_2)/d_2)/63.0;
        
        if( err < eps) {
            if(y_prev*y_22 < 0 && y_21 > 0) {
                h = h/2;
            }
            else {
                 fprintf(Out_t, "%5.11lf, ", t+h);
                 fprintf(Out_x, "%5.11lf, ", y_11);
                 fprintf(Out_y, "%5.11lf, ", y_12);
                 fprintf(Out_t, "%5.11lf, ", t + 2*h);
                 fprintf(Out_x, "%5.11lf, ", y_21);
                 fprintf(Out_y, "%5.11lf, ", y_22);
                y = y_22;
                x = y_21;
                t = t + 2*h;
                h = h*min(facmax, max(facmin, fac*pow(eps/err, 1.0/7.0 )));
                y_prev = y_22;
            }
        }
        else {
            h = h*min(1.0, max(facmin, fac*pow(eps/err, 1.0/7.0)) );
        }
        if( fabs(y_22) < eps && y_21 > 0) flag = false;
    }
    
    fclose(Out_t);
    fclose(Out_x);
    fclose(Out_y);
    
    return t;
}

void runge_numbers(double x, double t) {
    double y, q_11, q_12, q_21, q_22, q_31, q_32, tau;
    tau = 0;
    int i;
    y = 0;
    for(i = 3; i >= 0; i--) {
        q_11 = solution_1(x, 0, tau, t/(i+1), 1.e-8);
        q_12 = solution_2(x, 0, tau, t/(i+1), 1.e-8);
        q_21 = solution_1(x, 0, tau, t/(i+1), 1.e-10);
        q_22 = solution_2(x, 0, tau, t/(i+1), 1.e-10);
        q_31 = solution_1(x, 0, tau, t/(i+1), 1.e-12);
        q_32 = solution_2(x, 0,tau, t/(i+1), 1.e-12);
        printf("R_1 = %lf, R_2 = %lf, at i = %d\n",fabs((q_11 - q_21)/(q_21 - q_31)),fabs((q_12 - q_22)/(q_22 - q_32)), i);
    }
}

double solution_2(double x, double y, double t, double L, double eps) {
    double h, y_11, y_12, y_21, y_22, w_1, w_2, d_1, d_2;//y_prev;
    bool flag = true; // flag
    double err;
    
    h = 0.3;  // first step
    y_11 = 0;
    y_12 = 0;
    y_21 = 0;
    y_22 = 0;
    
    
    while(flag) {
        if(t + 2*h > L) h = 0.5*(L - t);
        y_11 = runge_1(x, y, t, h);
        y_12 = runge_2(x, y, t, h);
        y_21 = runge_1(y_11, y_12, t + h, h);
        y_22 = runge_2(y_11, y_12, t + h, h);
        w_1 = runge_1(x, y, t, 2*h);
        w_2 = runge_2(x, y, t, 2*h);
        d_1 = fabs(y_21 + (y_21 - w_1)/63);
        d_2 = fabs(y_22 + (y_22 - w_2)/63);
        err = max(fabs(y_21 - w_1)/d_1, fabs(y_22 - w_2)/d_2)/63.0;
        
        if( err < eps) {
            y = y_22;
            x = y_21;
            t = t + 2*h;
            h = h*min(facmax, max(facmin, fac*pow(eps/err, 1.0/7.0 )));
        }
        else {
            h = h*min(1.0, max(facmin, fac*pow(eps/err, 1.0/7.0)) );
        }
        if(fabs(L - t) < eps) flag = false;
    }
    return y;
}

double solution_1(double x, double y, double t, double L, double eps) {
    double h, y_11, y_12, y_21, y_22, w_1, w_2, d_1, d_2, delta;//y_prev;
    bool flag = true; // flag
    double err;
    
    h = 0.3;  // first step
    y_11 = 0;
    y_12 = 0;
    y_21 = 0;
    y_22 = 0;
    delta = 0;
    
    
    while(flag) {
        if(t + 2*h > L) h = 0.5*(L - t);
        y_11 = runge_1(x, y, t, h);
        y_12 = runge_2(x, y, t, h);
        y_21 = runge_1(y_11, y_12, t + h, h);
        y_22 = runge_2(y_11, y_12, t + h, h);
        w_1 = runge_1(x, y, t, 2*h);
        w_2 = runge_2(x, y, t, 2*h);
        d_1 = fabs(y_21 + (y_21 - w_1)/63);
        d_2 = fabs(y_22 + (y_22 - w_2)/63);
        err = max(fabs(y_21 - w_1)/d_1, fabs(y_22 - w_2)/d_2)/63.0;
        
        if( err < eps) {
            y = y_22;
            x = y_21;
            t = t + 2*h;
            delta = err + delta*e_lambdamax(x, y, h);
            h = h*min(facmax, max(facmin, fac*pow(eps/err, 1.0/7.0 )));
        }
        else {
            h = h*min(1.0, max(facmin, fac*pow(eps/err, 1.0/7.0)) );
        }
        if(fabs(L - t) < eps) flag = false;
    }
    printf("ERR_1 = %lf\n", delta);
    return x;
}

double e_lambdamax(double x, double y, double h) {
    return exp(max((1-3*alpha*x*x),(1 - 3*y*y))*h);
}
