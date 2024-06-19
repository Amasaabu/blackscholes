// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
// 
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice 
// Hall, John C. Hull,
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <thread>
#include <mutex>
#include <vector>

#include <time.h>
#include <chrono>
#include <iostream>

#include "thread_Pool.h"

//Precision to use for calculations
#define fptype float

#define NUM_RUNS 10000

typedef struct OptionData_ {
    fptype s;          // spot price
    fptype strike;     // strike price
    fptype r;          // risk-free interest rate
    fptype divq;       // dividend rate
    fptype v;          // volatility
    fptype t;          // time to maturity or option expiration in years 
    //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)  
    char OptionType;   // Option type.  "P"=PUT, "C"=CALL
    fptype divs;       // dividend vals (not used in this test)
    fptype DGrefval;   // DerivaGem Reference Value
} OptionData;

OptionData* data;
fptype* prices;
int numOptions;

int* otype;
fptype* sptprice;
fptype* strike;
fptype* rate;
fptype* volatility;
fptype* otime;
int numError = 0;

// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

fptype CNDF(fptype InputX)
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0) {
        InputX = -InputX;
        sign = 1;
    }
    else
        sign = 0;

    xInput = InputX;

    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;

    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal = xLocal_1 * xNPrimeofX;
    xLocal = 1.0 - xLocal;

    OutputX = xLocal;

    if (sign) {
        OutputX = 1.0 - OutputX;
    }

    return OutputX;
}

fptype BlkSchlsEqEuroNoDiv(fptype sptprice,
    fptype strike, fptype rate, fptype volatility,
    fptype time, int otype, float timet)
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1;
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    xSqrtTime = sqrt(xTime);

    logValues = log(sptprice / strike);

    xLogTerm = logValues;


    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 - xDen;

    d1 = xD1;
    d2 = xD2;

    NofXd1 = CNDF(d1);
    NofXd2 = CNDF(d2);

    FutureValueX = strike * (exp(-(rate) * (time)));
    if (otype == 0) {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    }
    else {
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}
fptype price;
std::condition_variable cv;
std::atomic<int> threads_finished_counter(0);
std::mutex waiter_mutex;

void dowrk(int start_option, int end_option) {
    for (int i = start_option; i < end_option; i++) {
        for (int j = 0; j < NUM_RUNS; j++) {
            price = BlkSchlsEqEuroNoDiv(sptprice[i], strike[i],
                rate[i], volatility[i], otime[i],
                otype[i], 0);
            prices[i] = price;
        }
    }
	//signal main that this thread is done
	{
		std::unique_lock<std::mutex> lock(waiter_mutex);
		threads_finished_counter++;
		cv.notify_one();
	}
};

int main(int argc, char** argv)
{
    FILE* file;
    int i;
    int loopnum;
    fptype* buffer;
    int* buffer2;
    int rv;

    if (argc != 3)
    {
        printf("Usage:\n\t%s <inputFile> <outputFile>\n", argv[0]);
        exit(1);
    }
    char* inputFile = argv[1];
    char* outputFile = argv[2];

    //Read input data from file
    file = fopen(inputFile, "r");
    if (file == NULL) {
        printf("ERROR: Unable to open file %s.\n", inputFile);
        exit(1);
    }
    rv = fscanf(file, "%i", &numOptions);
    if (rv != 1) {
        printf("ERROR: Unable to read from file %s.\n", inputFile);
        fclose(file);
        exit(1);
    }

    // alloc spaces for the option data
    data = (OptionData*)malloc(numOptions * sizeof(OptionData));
    prices = (fptype*)malloc(numOptions * sizeof(fptype));
    for (loopnum = 0; loopnum < numOptions; ++loopnum)
    {
        rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
        if (rv != 9) {
            printf("ERROR: Unable to read from file %s.\n", inputFile);
            fclose(file);
            exit(1);
        }
    }
    rv = fclose(file);
    if (rv != 0) {
        printf("ERROR: Unable to close file %s.\n", inputFile);
        exit(1);
    }

    printf("Num of Options: %d\n", numOptions);
    printf("Num of Runs: %d\n", NUM_RUNS);

#define PAD 256
#define LINESIZE 64

    buffer = (fptype*)malloc(5 * numOptions * sizeof(fptype) + PAD);
    sptprice = (fptype*)(((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = (int*)malloc(numOptions * sizeof(fptype) + PAD);
    otype = (int*)(((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));
    
    for (i = 0; i < numOptions; i++) {
        otype[i] = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i] = data[i].s;
        strike[i] = data[i].strike;
        rate[i] = data[i].r;
        volatility[i] = data[i].v;
        otime[i] = data[i].t;
    }


    printf("Size of data: %d\n", numOptions * (sizeof(OptionData) + sizeof(int)));


   //measuring time performanmce
    auto start = std::chrono::steady_clock::now();
    const int thread_count = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    int work_per_thread = numOptions / thread_count;

    //Thread_Pool thread_pool;
	for (int t = 0; t < thread_count; t++) {
		int start_option = t * work_per_thread;
		int end_option = (t == thread_count - 1) ? numOptions : start_option + work_per_thread;
        //thread_pool.submit([start_option, end_option]() {dowrk(start_option, end_option);});
		threads.push_back(std::thread(dowrk, start_option, end_option));
	}
    {
		std::cout << "Waiting for threads to finish...." << std::endl;
        std::unique_lock<std::mutex> lck(waiter_mutex);
        cv.wait(lck, [&]() {return threads_finished_counter==thread_count;});
    }

    for (auto& th : threads) {
        th.join();
    }
 
    //single threaded version
    //fptype price;
    //for (i = 0; i < numOptions; i++) {

    //    /* Calling main function to calculate option value based on
    //     * Black & Scholes's equation.
    //     */
    //    for (int j = 0; j < NUM_RUNS; j++) {
    //        price = BlkSchlsEqEuroNoDiv(sptprice[i], strike[i],
    //            rate[i], volatility[i], otime[i],
    //            otype[i], 0);
    //        prices[i] = price;
    //    }
    //}
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Operation completed in (milli): " << elapsed.count() << std::endl;
    std::cout << "writing to file now..." << std::endl;
    //Write prices to output file
    file = fopen(outputFile, "w");
    if (file == NULL) {
        printf("ERROR: Unable to open file %s.\n", outputFile);
        exit(1);
    }
    rv = fprintf(file, "%i\n", numOptions);
    if (rv < 0) {
        printf("ERROR: Unable to write to file %s.\n", outputFile);
        fclose(file);
        exit(1);
    }
    for (i = 0; i < numOptions; i++) {
        rv = fprintf(file, "%.18f\n", prices[i]);
        if (rv < 0) {
            printf("ERROR: Unable to write to file %s.\n", outputFile);
            fclose(file);
            exit(1);
        }
    }
    rv = fclose(file);
    if (rv != 0) {
        printf("ERROR: Unable to close file %s.\n", outputFile);
        exit(1);
    }

    free(data);
    free(prices);

    return 0;
}

