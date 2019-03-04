#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <tuple>
#include <cmath>

using std::ifstream;
using std::istream;
using std::ostream;
using std::getline;
using std::string;
using std::vector;
using std::tuple;

class windowTI {
public:
    windowTI() {}
    windowTI(istream& is);
    struct frameTI {
        unsigned long long timestep;
        double bond1;
        double avgbond1;
        double elect1;
        double avgelect1;
        double vdw1;
        double avgvdw1;
        double bond2;
        double avgbond2;
        double elect2;
        double avgelect2;
        double vdw2;
        double avgvdw2;
    };
    bool scanTILine(const char* str, frameTI& frame);
    void debugDump(ostream& os) const;
    ostream& compute(ostream& os) const;
    tuple<double, double, double, double, double, double> compute() const;
    friend void computeFreeEnergy(const vector<windowTI>& windows, ostream& os);
private:
    double mLambda;
    double mScalingBondP1;
    double mScalingBondP2;
    double mScalingVDWP1;
    double mScalingVDWP2;
    double mScalingELECP1;
    double mScalingELECP2;
    double mTemperature;
    vector<frameTI> mEquilFrames;
    vector<frameTI> mRunningFrames;
    unsigned long long mEquilSteps;
    unsigned long long mRunningSteps;
    unsigned long long mEquilNumLines;
    unsigned long long mRunningNumLines;
    unsigned long long mOutputFrequency;
};

windowTI::windowTI(istream& is) {
    string line;
    string lambda_info{"#NEW TI WINDOW: LAMBDA "};
    string partition1_info{"#PARTITION 1 SCALING: "};
    string partition2_info{"#PARTITION 2 SCALING: "};
    string temperature_info{"#CONSTANT TEMPERATURE: "};
    string equil_info{" STEPS OF EQUILIBRATION AT LAMBDA "};
    bool lambda_found = false;
    bool partition1_found = false;
    bool partition2_found = false;
    bool temperature_found = false;
    bool end_of_equilibrium = false;
    mEquilSteps = 0;
    mRunningSteps = 0;
    mEquilNumLines = 0;
    mRunningNumLines = 0;
    mOutputFrequency = 0;
    while (getline(is, line)) {
        if (!lambda_found) {
            string::size_type comp = line.find(lambda_info);
            if (comp != string::npos) {
                mLambda = std::stod(line.substr(comp + lambda_info.length(), line.length() - (comp + lambda_info.length())));
                lambda_found = true;
            }
        }
        if (!partition1_found) {
            string::size_type comp = line.find(partition1_info);
            if (comp != string::npos) {
                string bve_info = line.substr(comp + partition1_info.length(), line.length() - (comp + partition1_info.length()));
                std::sscanf(bve_info.c_str(), "BOND %lf VDW %lf ELEC %lf", &mScalingBondP1, &mScalingVDWP1, &mScalingELECP1);
                partition1_found = true;
            }
        }
        if (!partition2_found) {
            string::size_type comp = line.find(partition2_info);
            if (comp != string::npos) {
                string bve_info = line.substr(comp + partition2_info.length(), line.length() - (comp + partition2_info.length()));
                std::sscanf(bve_info.c_str(), "BOND %lf VDW %lf ELEC %lf", &mScalingBondP2, &mScalingVDWP2, &mScalingELECP2);
                partition2_found = true;
            }
        }
        if (!temperature_found) {
            string::size_type comp = line.find(temperature_info);
            if (comp != string::npos) {
                string t_K = line.substr(comp + temperature_info.length(), line.length() - (comp + temperature_info.length()));
                std::sscanf(t_K.c_str(), "%lf K", &mTemperature);
                temperature_found = true;
            }
        }
        if (!end_of_equilibrium) {
            frameTI tmp;
            if (scanTILine(line.c_str(), tmp)) {
                mEquilFrames.push_back(tmp);
                ++mEquilNumLines;
            } else {
                string::size_type comp = line.find(equil_info);
                if (comp != string::npos) {
                    double lambda_equil_line;
                    string scan_string = "#%llu" + equil_info + "%lf COMPLETED";
                    std::sscanf(line.c_str(), scan_string.c_str(), &mEquilSteps, &lambda_equil_line);
                    end_of_equilibrium = true;
                }
            }
        } else {
            frameTI tmp;
            if (scanTILine(line.c_str(), tmp)) {
                mRunningFrames.push_back(tmp);
                ++mRunningNumLines;
            }
        }
    }
    mRunningSteps = mRunningFrames.back().timestep;
    auto it = mRunningFrames.end() - 2;
    mOutputFrequency = mRunningSteps - (*it).timestep;
}

bool windowTI::scanTILine(const char* str, frameTI& frame) {
    if (std::sscanf(str, "TI: %llu %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &(frame.timestep), &(frame.bond1), &(frame.avgbond1), &(frame.elect1), &(frame.avgelect1), &(frame.vdw1), &(frame.avgvdw1), &(frame.bond2), &(frame.avgbond2), &(frame.elect2), &(frame.avgelect2), &(frame.vdw2), &(frame.avgvdw2)) == 13) {
        return true;
    } else {
        return false;
    }
}

void windowTI::debugDump(ostream& os) const {
    os << "Lambda = " << mLambda << '\n'
       << "P1 Scaling: Bond = " << mScalingBondP1 << " ; Vdw = " << mScalingVDWP1 << " ; Elect = " << mScalingELECP1 << '\n'
       << "P2 Scaling: Bond = " << mScalingBondP2 << " ; Vdw = " << mScalingVDWP2 << " ; Elect = " << mScalingELECP2 << '\n'
       << "Temperature = " << mTemperature << " K\n"
       << "Equilibration steps: " << mEquilSteps << " ; lines = " << mEquilNumLines << '\n'
       << "Running steps: " << mRunningSteps << " ; lines = " << mRunningNumLines << '\n'
       << "OutputFrequency: " << mOutputFrequency << '\n';
}

ostream& windowTI::compute(ostream& os) const {
    std::ios_base::fmtflags f(os.flags());
    int output_width = 12;
    os << std::fixed << std::setprecision(5);
    os << mLambda << " ";
    auto [mean_bond, var_bond, mean_vdw, var_vdw, mean_elect, var_elect] = compute();
    os << std::setw(output_width) << mean_bond << " " << std::setw(output_width) << var_bond << " "
       << std::setw(output_width) << mean_vdw << " " << std::setw(output_width) << var_vdw << " "
       << std::setw(output_width) << mean_elect << " " << std::setw(output_width) << var_elect << " ";
    os.flags(f);
    return os;
}

tuple<double, double, double, double, double, double> windowTI::compute() const {
    double dUdl_bond = 0, dUdl_vdw = 0, dUdl_elect = 0;
    double dUdl_bond2 = 0, dUdl_vdw2 = 0, dUdl_elect2 = 0;
    for (size_t i = 0; i < mRunningFrames.size(); ++i) {
        dUdl_bond += mRunningFrames[i].bond1 - mRunningFrames[i].bond2;
        dUdl_vdw += mRunningFrames[i].vdw1 - mRunningFrames[i].vdw2;
        dUdl_elect += mRunningFrames[i].elect1 - mRunningFrames[i].elect2;
        dUdl_bond2 += (mRunningFrames[i].bond1 - mRunningFrames[i].bond2) * (mRunningFrames[i].bond1 - mRunningFrames[i].bond2);
        dUdl_vdw2 += (mRunningFrames[i].vdw1 - mRunningFrames[i].vdw2) * (mRunningFrames[i].vdw1 - mRunningFrames[i].vdw2);
        dUdl_elect2 += (mRunningFrames[i].elect1 - mRunningFrames[i].elect2) * (mRunningFrames[i].elect1 - mRunningFrames[i].elect2);
    }
    const double mean_bond = dUdl_bond / mRunningNumLines;
    const double mean_vdw = dUdl_vdw / mRunningNumLines;
    const double mean_elect = dUdl_elect / mRunningNumLines;
    const double var_bond = (dUdl_bond2 - mean_bond * mean_bond * mRunningNumLines) / (mRunningNumLines - 1);
    const double var_vdw = (dUdl_vdw2 - mean_vdw * mean_vdw * mRunningNumLines) / (mRunningNumLines - 1);
    const double var_elect = (dUdl_elect2 - mean_elect * mean_elect * mRunningNumLines) / (mRunningNumLines - 1);
    return std::make_tuple(mean_bond, var_bond, mean_vdw, var_vdw, mean_elect, var_elect);
}

void computeFreeEnergy(const vector<windowTI>& windows, ostream& os) {
    int output_width = 12;
    double dF = 0, sigma_dF = 0;
    std::ios_base::fmtflags f(os.flags());
    os << std::fixed << std::setprecision(5);
    os << windows[0].mLambda;
    os << std::setw(output_width) << dF
        << std::setw(output_width) << sigma_dF
        << '\n';
    for (size_t i = 1; i < windows.size(); ++i) {
        auto [mean_bond, var_bond, mean_vdw, var_vdw, mean_elect, var_elect] = windows[i-1].compute();
        auto [mean_bond_next, var_bond_next, mean_vdw_next, var_vdw_next, mean_elect_next, var_elect_next] = windows[i].compute();
        dF += 0.5 * (mean_bond_next + mean_bond + mean_vdw_next + mean_vdw + mean_elect_next + mean_elect) * (windows[i].mLambda - windows[i-1].mLambda);
        sigma_dF += 0.25 * (var_bond_next + var_bond + var_vdw_next + var_vdw + var_elect_next + var_elect) * (windows[i].mLambda - windows[i-1].mLambda) * (windows[i].mLambda - windows[i-1].mLambda);
        os << windows[i].mLambda;
        os << std::setw(output_width) << dF
           << std::setw(output_width) << std::sqrt(sigma_dF)
           << '\n';
    }
    os.flags(f);
}

int main(int argc, char* argv[]) {
    ifstream ifsTIOutput(argv[1]);
    if (!ifsTIOutput.is_open()) {
        std::cerr << "Error opening file!" << '\n';
        return 1;
    }
    std::stringstream data;
    string line;
    bool first_time = true;
    vector<windowTI> TIresults;
    while (std::getline(ifsTIOutput, line)) {
        if (line.find(std::string{"#NEW TI WINDOW:"}) != string::npos) {
            if (!first_time) {
                TIresults.push_back(windowTI(data));
            }
            first_time = false;
            data.clear();
            data << line << '\n';
        } else if (line.find(std::string{"#TITITLE:"}) != string::npos) {
            continue;
        } else {
            data << line << '\n';
        }
    }
    bool output_mean_variance = true;
    if (output_mean_variance) {
        std::cout << "     #λ   Mean(bond)    Var(bond)    Mean(vdw)     Var(vdw)  Mean(elect)   Var(elect)\n";
        for (const auto& i : TIresults) {
            i.compute(std::cout) << std::endl;
        }
    }
    TIresults.push_back(windowTI(data));
    std::cout << "\n     #λ       dF(λ)     σ_dF(λ)" << std::endl;
    computeFreeEnergy(TIresults, std::cout);
    return 0;
}
