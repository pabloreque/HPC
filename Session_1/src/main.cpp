// Numerically demonstrates the order of accuracy of two finite difference
// formulas for the first derivative:
//   forward: (f(x+dx) - f(x)) / dx            -> O(dx)   (1st order)
//   central: (f(x+dx) - f(x-dx)) / (2*dx)     -> O(dx^2) (2nd order)
//
// Uses f(x) = sin(x), with exact derivative f'(x) = cos(x), and measures
// the absolute error for a sequence of step sizes dx that is halved on
// every iteration. The observed order is estimated as:
//   p = log2( error(dx) / error(dx/2) )
// because halving dx at each step means an error ~ C*dx^p implies
// error(dx)/error(dx/2) = 2^p.

#include <H5Cpp.h>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>

double forwardDiff(const std::function<double(double)> &f, double x, double dx) {
    return (f(x + dx) - f(x)) / dx;
}

double centralDiff(const std::function<double(double)> &f, double x, double dx) {
    return (f(x + dx) - f(x - dx)) / (2.0 * dx);
}

int main() {
    const double x0 = 1.0;
    auto f = [](double x) { return std::sin(x); };
    const double exact = std::cos(x0);

    const int numSteps = 15;
    const double dx0 = 0.1;

    std::vector<double> dxs(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        dxs[i] = dx0 / std::pow(2.0, i);
    }

    std::vector<double> errForward(numSteps), errCentral(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        errForward[i] = std::fabs(forwardDiff(f, x0, dxs[i]) - exact);
        errCentral[i] = std::fabs(centralDiff(f, x0, dxs[i]) - exact);
    }

    std::printf("%-12s %-14s %-10s %-14s %-10s\n", "dx", "err_forward", "p_forward", "err_central", "p_central");
    for (int i = 0; i < numSteps; ++i) {
        std::printf("%-12.6e %-14.6e", dxs[i], errForward[i]);
        if (i == 0) {
            std::printf(" %-10s", "-");
        } else {
            double p = std::log2(errForward[i - 1] / errForward[i]);
            std::printf(" %-10.4f", p);
        }
        std::printf(" %-14.6e", errCentral[i]);
        if (i == 0) {
            std::printf(" %-10s", "-");
        } else {
            double p = std::log2(errCentral[i - 1] / errCentral[i]);
            std::printf(" %-10.4f", p);
        }
        std::printf("\n");
    }

    const std::string outFile = "results.h5";
    H5::H5File file(outFile, H5F_ACC_TRUNC);
    hsize_t dims[1] = {static_cast<hsize_t>(numSteps)};
    H5::DataSpace dataspace(1, dims);

    H5::DataSet dxDataset = file.createDataSet("dx", H5::PredType::NATIVE_DOUBLE, dataspace);
    dxDataset.write(dxs.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet errForwardDataset = file.createDataSet("err_forward", H5::PredType::NATIVE_DOUBLE, dataspace);
    errForwardDataset.write(errForward.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet errCentralDataset = file.createDataSet("err_central", H5::PredType::NATIVE_DOUBLE, dataspace);
    errCentralDataset.write(errCentral.data(), H5::PredType::NATIVE_DOUBLE);

    std::cout << "\nResults exported to " << outFile << "\n";
    std::cout << "p_forward should converge to ~1.0 and p_central to ~2.0\n";
    std::cout << "as dx decreases (until round-off error starts to dominate).\n";

    return 0;
}
