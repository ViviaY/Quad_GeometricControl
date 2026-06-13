# quadrotor-safe-synthesis

Code for the submission Reach-Avoid Control Synthesis for a Quadrotor UAV with Formal Safety Guarantees

Written on MATLAB R2023a


Dependencies:
1. MATLAB Global Optimization Toolbox (for simulated annealing)
2. CasADi v3.7.1 (for nonlinear optimization and QP solving): https://web.casadi.org


This project consists of three main contollers:
    1. The proposed controller Geometric control:
        main_GC.m: regenerating all models and control gains with different parameter settings;
        fast_main.m: using the paper's parameters and models to reproduce the results.
        cpu_check.m：testing the runtimes for different paramters in RRT, the results are shown in Table 1.
    
    2. Baseline 1 control barrier function[1]: main_CBFQP.m is used to track the reference trajecotry.

    3. Baseline 2 Nonlinear model predictive control on SE(3)[2,3]: main_NMPC.m is used to track the reference trajecotry;


SafetyReeachingChecking.m: used to check the reference trajecotry safety and reach margin, and the safety and reach rates for the controller's tracking results, the results are shown in Table 2. 

plotting_python: plotting the Figures 2-4 of the paper. 

Note that all results reported in the paper are saved as .mat files in the controllers’ subfolders named "results_submit".


[1] V. Freire and X. Xu, “Flatness-based quadcopter trajectory planning and tracking with continuous-time safety guarantees,” IEEE Transactions on Control Systems Technology, vol. 31, no. 6, pp. 2319–2334, 2023.
[2] M. Elhesasy, T. N. Dief, M. Atallah, M. Okasha, M. M. Kamra, S. Yoshida, and M. A. Rushdi, "Non-linear model predictive control using casadi package for trajectory tracking of quadrotor," Energies, vol. 16, no. 5, p. 2143, 2023.
[3] J. C. Pereira, V. J. Leite, and G. V. Raffo, "Nonlinear model predictive control on SE(3) for quadrotor aggressive maneuvers," Journal of Intelligent \& Robotic Systems, vol. 101, no. 3, p. 62, 2021.
