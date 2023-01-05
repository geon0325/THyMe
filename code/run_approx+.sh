dataset=email-Enron-full
delta=86400000
T=0.001
S=100

g++ -O3 -std=c++11 main_approx+.cpp -o run_approx+;
./run_approx+ $dataset $delta $T $S;
rm run_approx+;
