g++ -O3 -std=c++11 main_thyme.cpp -o run_thyme;

dataset=email-Enron-full
delta=86400000
./run_thyme $dataset $delta;
rm run_thyme;