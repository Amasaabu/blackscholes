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
#include "Mpi_Lib.h"
#include <time.h>
#include <chrono>
#include <iostream>
#include "mpi.h"

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

OptionData* data_;
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


//void dowork(int start_option, int end_option) {
//    for (int i = start_option; i < end_option; i++) {
//    std::cout << "Running work/..." << " Start option " << start_option << " End Option" << end_option << std::endl;
//    for (int j = 0; j < NUM_RUNS; j++) {
//        std::cout << NUM_RUNS << std::endl;
//        price = BlkSchlsEqEuroNoDiv(sptprice[i], strike[i],
//            rate[i], volatility[i], otime[i],
//            otype[i], 0);
//        prices[i] = price;
//        std::cout << "Price: " << price << std::endl;
//    }
//}
//}
//void dowrk(int numbs) {
//    for (int i = 0; i < numbs; i++) {
//
//        /* Calling main function to calculate option value based on
//         * Black & Scholes's equation.
//         */
//        for (int j = 0; j < NUM_RUNS; j++) {
//            price = BlkSchlsEqEuroNoDiv(sptprice[i], strike[i],
//                rate[i], volatility[i], otime[i],
//                otype[i], 0);
//            prices[i] = price;
//        }
//    }
//};

int main(int argc, char** argv)
{

    Mpi_Lib<OptionData, int, float> mpi_lib(argc, argv);
    int world_rank = mpi_lib.get_world_rank();
    int i;
    int loopnum;
    fptype* buffer;
    int* buffer2;
    char* outputFile;
    FILE* file;
    int rv;

    if (world_rank == 0) {
        if (argc != 3)
        {
            printf("Usage:\n\t%s <inputFile> <outputFile>\n", argv[0]);
            exit(1);
        }
        char* inputFile = argv[1];
        outputFile = argv[2];

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
        data_ = (OptionData*)malloc(numOptions * sizeof(OptionData));
        prices = (fptype*)malloc(numOptions * sizeof(fptype));
        for (loopnum = 0; loopnum < numOptions; ++loopnum)
        {
            rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data_[loopnum].s, &data_[loopnum].strike, &data_[loopnum].r, &data_[loopnum].divq, &data_[loopnum].v, &data_[loopnum].t, &data_[loopnum].OptionType, &data_[loopnum].divs, &data_[loopnum].DGrefval);
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
            otype[i] = (data_[i].OptionType == 'P') ? 1 : 0;
            sptprice[i] = data_[i].s;
            strike[i] = data_[i].strike;
            rate[i] = data_[i].r;
            volatility[i] = data_[i].v;
            otime[i] = data_[i].t;
        }


        printf("Size of datas: %d\n", numOptions * (sizeof(OptionData) + sizeof(int)));
    }
   // mpi_lib.init(numOptions, sizeof(OptionData) + sizeof(int));
    mpi_lib.init(numOptions, sizeof(OptionData));

    mpi_lib.broadcast(&numOptions, 1, 0);

   
   auto sendcounts = mpi_lib.get_sendcounts();
  //  int num_local_options = sendcounts[world_rank] / sizeof(OptionData);
    
    OptionData* local_data=new OptionData[sendcounts[world_rank]];
    float *local_results= new float [sendcounts[world_rank]];

    //ensure every processor has work to do
    if (sendcounts[world_rank] == 0) {
        cout << "Processor " << world_rank << " has no work to do" << endl;
        MPI_Finalize();
        return 0;
    }
  //  std::cout << "Processor 0 distributing load...." << std::endl;
    mpi_lib.scatterV(data_, local_data, [ & local_results, &local_data](int start, int end) {
        std::cout << "Start " << start << " End " << end << std::endl;
        for (int i = start; i < end; i++) {
            // Perform calculations for each option
            std::cout << "Running work/..." << " Start option " << start << " End Option" << end << std::endl;
            for (int j = 0; j < NUM_RUNS; j++) {
                local_results[i] = BlkSchlsEqEuroNoDiv(local_data[i].s, local_data[i].strike,
                    local_data[i].r, local_data[i].v, local_data[i].t,
                    (local_data[i].OptionType == 'P' ? 1 : 0), 0);
                //prices[i] = price;
            }
            std::cout << "Price: " << local_results[0] << std::endl;
        }
    });
    std::cout << "Processor 0 Finished distributing load...." << std::endl;

   
   // mpi_lib.barrier();
    mpi_lib.gather_v(local_results, prices);
    if (world_rank == 0) {
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

        free(data_);
        free(prices);
    }
    return 0;
}
