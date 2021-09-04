g++ -O3 -std=c++11 main_dp.cpp -o run_dp;

dataset=email-Enron-full
for delta in 86400000
do
	for random in 0 1
	do
		./run_dp $dataset $random $delta;
	done
done
rm run_dp;

