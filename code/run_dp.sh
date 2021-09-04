g++ -O3 -std=c++11 main_dp.cpp -o run_dp;

dataset=email-Enron-full
delta=86400000
./run_dp $dataset $delta;
rm run_dp

