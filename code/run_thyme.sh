g++ -O3 -std=c++11 main_thyme.cpp -o run_thyme;

dataset=email-Enron-full
for delta in 86400000
do
	for random in 0 1
	do
		./run_thyme $dataset $random $delta;
	done
done
rm run_thyme;

