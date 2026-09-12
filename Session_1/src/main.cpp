// Numerically demonstrates the order of accuracy of two finite difference
// formulas for the first derivative:
//   forward: (f(x+dx) - f(x)) / dx            -> O(dx)   (1st order)
//   central: (f(x+dx) - f(x-dx)) / (2*dx)     -> O(dx^2) (2nd order)
//
// Uses f(x) = sin(x), with exact derivative f'(x) = cos(x), and measures
// the absolute error for a sequence of step sizes dx that is halved on
// every iteration. The observed order is estimated as:
//   p = log2( error(dx) / error(dx/2) )
// This is because the error is ~ C*dx^p, which implies
//   error(dx)/error(dx/2) = 2^p.

#include <H5Cpp.h>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

// Named alias so signatures read as "takes a scalar function" instead of the raw std::function type.
using ScalarFunction = std::function<double(double)>;

// Fixed-field bundle (like a Python dataclass) so the three parallel vectors travel
// together as one coherent result instead of three loose arguments.
struct ConvergenceResults {
    std::vector<double> stepSizes;
    std::vector<double> forwardErrors;
    std::vector<double> centralErrors;
};

double forwardDiff(const ScalarFunction &f, double x, double dx) {
    return (f(x + dx) - f(x)) / dx;
}

double centralDiff(const ScalarFunction &f, double x, double dx) {
    return (f(x + dx) - f(x - dx)) / (2.0 * dx);
}

std::vector<double> generateHalvingStepSizes(double initialStep, int numSteps) {
    std::vector<double> stepSizes(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        stepSizes[i] = initialStep / std::pow(2.0, i);
    }
    return stepSizes;
}

std::vector<double> computeAbsoluteErrors(const ScalarFunction &f, double evaluationPoint, double exactDerivative,
                                           const std::vector<double> &stepSizes,
                                           const std::function<double(const ScalarFunction &, double, double)> &diffMethod) {
    std::vector<double> errors(stepSizes.size());
    for (std::size_t i = 0; i < stepSizes.size(); ++i) {
        errors[i] = std::fabs(diffMethod(f, evaluationPoint, stepSizes[i]) - exactDerivative);
    }
    return errors;
}

ConvergenceResults runConvergenceStudy(const ScalarFunction &f, double evaluationPoint, double exactDerivative,
                                        double initialStep, int numSteps) {
    ConvergenceResults results;
    results.stepSizes = generateHalvingStepSizes(initialStep, numSteps);
    results.forwardErrors = computeAbsoluteErrors(f, evaluationPoint, exactDerivative, results.stepSizes, forwardDiff);
    results.centralErrors = computeAbsoluteErrors(f, evaluationPoint, exactDerivative, results.stepSizes, centralDiff);
    return results;
}

double computeObservedOrder(double previousError, double currentError) {
    return std::log2(previousError / currentError);
}

void printOrderColumn(const std::vector<double> &errors, std::size_t i) {
    if (i == 0) {
        std::printf(" %-10s", "-");
    } else {
        std::printf(" %-10.4f", computeObservedOrder(errors[i - 1], errors[i]));
    }
}

void printConvergenceTable(const ConvergenceResults &results) {
    std::printf("%-12s %-14s %-10s %-14s %-10s\n", "dx", "err_forward", "p_forward", "err_central", "p_central");

    for (std::size_t i = 0; i < results.stepSizes.size(); ++i) {
        std::printf("%-12.6e %-14.6e", results.stepSizes[i], results.forwardErrors[i]);
        printOrderColumn(results.forwardErrors, i);
        std::printf(" %-14.6e", results.centralErrors[i]);
        printOrderColumn(results.centralErrors, i);
        std::printf("\n");
    }
}

void writeDataset(H5::H5File &file, const std::string &name, const std::vector<double> &data,
                   const H5::DataSpace &dataspace) {
    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void exportResultsToHDF5(const std::string &outFile, const ConvergenceResults &results) {
    H5::H5File file(outFile, H5F_ACC_TRUNC);
    hsize_t dims[1] = {static_cast<hsize_t>(results.stepSizes.size())};
    H5::DataSpace dataspace(1, dims);

    writeDataset(file, "dx", results.stepSizes, dataspace);
    writeDataset(file, "err_forward", results.forwardErrors, dataspace);
    writeDataset(file, "err_central", results.centralErrors, dataspace);
}

void printExportSummary(const std::string &outFile) {
    std::cout << "\nResults exported to " << outFile << "\n";
    std::cout << "p_forward should converge to ~1.0 and p_central to ~2.0\n";
    std::cout << "as dx decreases (until round-off error starts to dominate).\n";
}

int main() {
    const double evaluationPoint = 1.0;
    const auto f = [](double x) { return std::sin(x); };
    const double exactDerivative = std::cos(evaluationPoint);
    const int numSteps = 15;
    const double dx0 = 0.1;
    const std::string outFile = "results.h5";

    const ConvergenceResults results = runConvergenceStudy(f, evaluationPoint, exactDerivative, dx0, numSteps);

    printConvergenceTable(results);
    exportResultsToHDF5(outFile, results);
    printExportSummary(outFile);

    return 0;
}
