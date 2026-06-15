#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>

// Adjust include path as needed
#include "mcfadbad.hpp"

using std::cout;
using std::endl;

using FADF = fadbad::F<double>;
using BADB = fadbad::B<double>;

static constexpr double PI =
    3.141592653589793238462643383279502884;

// ============================================================
// Utilities
// ============================================================

double NaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool is_close(double a, double b, double atol = 1e-10, double rtol = 1e-8)
{
    if (std::isnan(a) && std::isnan(b)) return true;
    if (std::isnan(a) || std::isnan(b)) return false;
    if (std::isinf(a) || std::isinf(b)) return a == b;
    const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    return std::fabs(a - b) <= atol + rtol * scale;
}

void print_check(double ad, double ref, double atol = 1e-10, double rtol = 1e-8)
{
    const bool ok = is_close(ad, ref, atol, rtol);
    const double err = (std::isnan(ad) || std::isnan(ref)) ? NaN() : std::fabs(ad - ref);

    cout << std::setw(16) << ad
         << std::setw(16) << ref
         << std::setw(16) << err
         << std::setw(10) << (ok ? "OK" : "FAIL");
}

std::vector<double> make_grid(double a, double b, double h)
{
    std::vector<double> g;
    for (double x = a; x <= b + 0.5*h; x += h) g.push_back(x);
    return g;
}

// ============================================================
// Analytical reference functions
// ============================================================

// ---- abs(x) ----
double abs_val(double x)
{
    return std::fabs(x);
}

double abs_der(double x)
{
    if (x < 0.0) return -1.0;
    if (x > 0.0) return  1.0;
    return NaN(); // non-differentiable at 0
}

// ---- erf(x) ----
double erf_val(double x)
{
    return std::erf(x);
}

double erf_der(double x)
{
    return 2.0 / std::sqrt(PI) * std::exp(-x*x);
}

// ---- xlog(x) = x log x ----
double xlog_val(double x)
{
    return x * std::log(x);
}

double xlog_der(double x)
{
    return std::log(x) + 1.0;
}

// ---- max(x,y) ----
double max_val(double x, double y)
{
    return std::max(x, y);
}

double dmax_dx(double x, double y)
{
    if (x > y) return 1.0;
    if (x < y) return 0.0;
    return NaN(); // tie: depends on your implementation
}

double dmax_dy(double x, double y)
{
    if (x > y) return 0.0;
    if (x < y) return 1.0;
    return NaN(); // tie: depends on your implementation
}

// ---- min(x,y) ----
double min_val(double x, double y)
{
    return std::min(x, y);
}

double dmin_dx(double x, double y)
{
    if (x < y) return 1.0;
    if (x > y) return 0.0;
    return NaN(); // tie: depends on your implementation
}

double dmin_dy(double x, double y)
{
    if (x < y) return 0.0;
    if (x > y) return 1.0;
    return NaN(); // tie: depends on your implementation
}

// ---- max(x,c), min(x,c) with c = 1.0 ----
double max1_val(double x)
{
    return std::max(x, 1.0);
}

double dmax1_dx(double x)
{
    if (x > 1.0) return 1.0;
    if (x < 1.0) return 0.0;
    return NaN(); // non-differentiable at x = 1
}

double min1_val(double x)
{
    return std::min(x, 1.0);
}

double dmin1_dx(double x)
{
    if (x < 1.0) return 1.0;
    if (x > 1.0) return 0.0;
    return NaN(); // non-differentiable at x = 1
}
// ============================================================
// Extra helper derivatives for shifted / scalar-kink functions
// ============================================================

double sign_shift(double x, double c)
{
    if (x < c) return -1.0;
    if (x > c) return  1.0;
    return NaN();
}

double dmax_scalar_dx(double x, double c)
{
    if (x > c) return 1.0;
    if (x < c) return 0.0;
    return NaN();
}

double dmin_scalar_dx(double x, double c)
{
    if (x < c) return 1.0;
    if (x > c) return 0.0;
    return NaN();
}

// ============================================================
// Composite unary test 1
//
// f(x) = erf(x) * xlog(x+2) + 0.5*abs(x-0.3)
// with xlog(z) = z log(z)
// domain: x > -2
// ============================================================

double comp1_val(double x)
{
    const double z = x + 2.0;
    return std::erf(x) * (z * std::log(z)) + 0.5 * std::fabs(x - 0.3);
}

double comp1_der(double x)
{
    const double z = x + 2.0;
    const double dz = std::log(z) + 1.0;
    const double der_erf = 2.0 / std::sqrt(PI) * std::exp(-x * x);
    return der_erf * (z * std::log(z))
         + std::erf(x) * dz
         + 0.5 * sign_shift(x, 0.3);
}

// ============================================================
// Composite unary test 2
//
// f(x) = max(x,1.0) * min(x+0.5,1.0) + x^2
// kinks at x = 1.0 and x = 0.5
// ============================================================

double comp2_val(double x)
{
    const double a = std::max(x, 1.0);
    const double b = std::min(x + 0.5, 1.0);
    return a * b + x * x;
}

double comp2_der(double x)
{
    const double a  = std::max(x, 1.0);
    const double b  = std::min(x + 0.5, 1.0);
    const double da = dmax_scalar_dx(x, 1.0);
    const double db = dmin_scalar_dx(x + 0.5, 1.0);  // derivative wrt x

    if (std::isnan(da) || std::isnan(db)) return NaN();
    return da * b + a * db + 2.0 * x;
}

// ============================================================
// Composite binary test 3
//
// f(x,y) = max(x,y)*(x+y) + min(x,1.0)*y^2 + erf(x-y)
// kinks at x=y and x=1
// ============================================================

double comp3_val(double x, double y)
{
    const double a = std::max(x, y);
    const double b = std::min(x, 1.0);
    return a * (x + y) + b * y * y + std::erf(x - y);
}

double comp3_dx(double x, double y)
{
    const double a   = std::max(x, y);
    const double da  = dmax_dx(x, y);
    const double bdx = dmin_scalar_dx(x, 1.0);
    const double der_erf = 2.0 / std::sqrt(PI) * std::exp(-(x - y) * (x - y));

    if (std::isnan(da) || std::isnan(bdx)) return NaN();
    return da * (x + y) + a + bdx * y * y + der_erf;
}

double comp3_dy(double x, double y)
{
    const double a   = std::max(x, y);
    const double da  = dmax_dy(x, y);
    const double b   = std::min(x, 1.0);
    const double der_erf = 2.0 / std::sqrt(PI) * std::exp(-(x - y) * (x - y));

    if (std::isnan(da)) return NaN();
    return da * (x + y) + a + 2.0 * b * y - der_erf;
}

// ============================================================
// Composite binary test 4
//
// f(x,y) = abs(x-0.5)*max(y,1.0) + xlog(y+2.5) + x*y
// with xlog(z) = z log(z)
// kinks at x = 0.5 and y = 1
// domain: y > -2.5
// ============================================================

double comp4_val(double x, double y)
{
    const double z = y + 2.5;
    return std::fabs(x - 0.5) * std::max(y, 1.0)
         + z * std::log(z)
         + x * y;
}

double comp4_dx(double x, double y)
{
    return sign_shift(x, 0.5) * std::max(y, 1.0) + y;
}

double comp4_dy(double x, double y)
{
    const double z = y + 2.5;
    const double dmax = dmax_scalar_dx(y, 1.0);
    if (std::isnan(dmax)) return NaN();

    return std::fabs(x - 0.5) * dmax
         + (std::log(z) + 1.0)
         + x;
}

// ============================================================
// Forward-mode tests
// ============================================================

template <typename FADFunc, typename ValFunc, typename DerFunc>
void test_unary_F_compare(
    const std::string& name,
    FADFunc fad_f,
    ValFunc ref_f,
    DerFunc ref_df,
    const std::vector<double>& grid)
{
    cout << "\n================================================================================================================\n";
    cout << "F<double> unary test: " << name << "\n";
    cout << "================================================================================================================\n";
    cout << std::setw(12) << "x"
         << std::setw(16) << "AD f"
         << std::setw(16) << "REF f"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check"
         << std::setw(16) << "AD df"
         << std::setw(16) << "REF df"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check" << "\n";

    for (double xv : grid) {
        FADF x = xv;
        x.diff(0, 1);
        FADF y = fad_f(x);

        cout << std::setw(12) << xv;
        print_check(y.val(), ref_f(xv));
        print_check(y[0],  ref_df(xv));
        cout << "\n";
    }
}

template <typename FADFunc, typename ValFunc, typename DxFunc, typename DyFunc>
void test_binary_F_compare(
    const std::string& name,
    FADFunc fad_f,
    ValFunc ref_f,
    DxFunc ref_dx,
    DyFunc ref_dy,
    const std::vector<double>& xgrid,
    const std::vector<double>& ygrid)
{
    cout << "\n============================================================================================================================================================\n";
    cout << "F<double> binary test: " << name << "\n";
    cout << "============================================================================================================================================================\n";
    cout << std::setw(10) << "x"
         << std::setw(10) << "y"
         << std::setw(16) << "AD f"
         << std::setw(16) << "REF f"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check"
         << std::setw(16) << "AD dfdx"
         << std::setw(16) << "REF dfdx"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check"
         << std::setw(16) << "AD dfdy"
         << std::setw(16) << "REF dfdy"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check" << "\n";

    for (double xv : xgrid) {
        for (double yv : ygrid) {
            FADF x = xv; x.diff(0, 2);
            FADF y = yv; y.diff(1, 2);
            FADF z = fad_f(x, y);

            cout << std::setw(10) << xv
                 << std::setw(10) << yv;
            print_check(z.val(), ref_f(xv, yv));
            print_check(z[0],    ref_dx(xv, yv));
            print_check(z[1],    ref_dy(xv, yv));
            cout << "\n";
        }
        cout << "\n";
    }
}

// ============================================================
// Reverse-mode tests
// ============================================================

template <typename BADFunc, typename ValFunc, typename DerFunc>
void test_unary_B_compare(
    const std::string& name,
    BADFunc bad_f,
    ValFunc ref_f,
    DerFunc ref_df,
    const std::vector<double>& grid)
{
    cout << "\n================================================================================================================\n";
    cout << "B<double> unary test: " << name << "\n";
    cout << "================================================================================================================\n";
    cout << std::setw(12) << "x"
         << std::setw(16) << "AD f"
         << std::setw(16) << "REF f"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check"
         << std::setw(16) << "AD df"
         << std::setw(16) << "REF df"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check" << "\n";

    for (double xv : grid) {
        BADB x = xv;
        BADB y = bad_f(x);
        y.diff(0, 1);

        cout << std::setw(12) << xv;
        print_check(y.val(), ref_f(xv));
        print_check(x.d(0),  ref_df(xv));
        cout << "\n";
    }
}

template <typename BADFunc, typename ValFunc, typename DxFunc, typename DyFunc>
void test_binary_B_compare(
    const std::string& name,
    BADFunc bad_f,
    ValFunc ref_f,
    DxFunc ref_dx,
    DyFunc ref_dy,
    const std::vector<double>& xgrid,
    const std::vector<double>& ygrid)
{
    cout << "\n============================================================================================================================================================\n";
    cout << "B<double> binary test: " << name << "\n";
    cout << "============================================================================================================================================================\n";
    cout << std::setw(10) << "x"
         << std::setw(10) << "y"
         << std::setw(16) << "AD f"
         << std::setw(16) << "REF f"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check"
         << std::setw(16) << "AD dfdx"
         << std::setw(16) << "REF dfdx"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check"
         << std::setw(16) << "AD dfdy"
         << std::setw(16) << "REF dfdy"
         << std::setw(16) << "|err|"
         << std::setw(10) << "check" << "\n";

    for (double xv : xgrid) {
        for (double yv : ygrid) {
            BADB x = xv;
            BADB y = yv;
            BADB z = bad_f(x, y);
            z.diff(0, 1);

            cout << std::setw(10) << xv
                 << std::setw(10) << yv;
            print_check(z.val(), ref_f(xv, yv));
            print_check(x.d(0),  ref_dx(xv, yv));
            print_check(y.d(0),  ref_dy(xv, yv));
            cout << "\n";
        }
        cout << "\n";
    }
}

// ============================================================
// Main
// ============================================================

int main()
{
    cout << std::fixed << std::setprecision(8);

    const auto grid_sym = make_grid(-2.0, 2.0, 0.5);
    const auto grid_pos = make_grid(0.25, 3.00, 0.25);
    const std::vector<double> grid2 = { -1.0, -0.25, 0.0, 0.25, 1.0 };

    const std::vector<double> grid_comp_u1 = { -1.5, -0.7, 0.0, 0.8, 1.6 };
    const std::vector<double> grid_comp_u2 = { -1.2, -0.2, 0.2, 0.8, 1.4 };

    const std::vector<double> grid_comp_x  = { -1.0, 0.4, 1.3 };
    const std::vector<double> grid_comp_y  = { -0.5, 0.2, 1.7 };

    const std::vector<double> grid_comp_x2 = { -1.0, 0.0, 1.4 };
    const std::vector<double> grid_comp_y2 = { -1.2, 0.3, 1.6 };

    // -----------------------------
    // F<double>
    // -----------------------------
    test_unary_F_compare(
        "fabs(x)",
        [](const FADF& x) { return fadbad::fabs(x); },
        abs_val, abs_der, grid_sym
    );

    test_unary_F_compare(
        "erf(x)",
        [](const FADF& x) { return fadbad::erf(x); },
        erf_val, erf_der, grid_sym
    );

    test_unary_F_compare(
        "xlog(x)  [assumed x*log(x)]",
        [](const FADF& x) { return fadbad::xlog(x); },
        xlog_val, xlog_der, grid_pos
    );

    test_binary_F_compare(
        "max(x,y)",
        [](const FADF& x, const FADF& y) { return fadbad::max(x, y); },
        max_val, dmax_dx, dmax_dy, grid2, grid2
    );

    test_binary_F_compare(
        "min(x,y)",
        [](const FADF& x, const FADF& y) { return fadbad::min(x, y); },
        min_val, dmin_dx, dmin_dy, grid2, grid2
    );

    test_unary_F_compare(
        "max(x,1.0)",
        [](const FADF& x) { return fadbad::max(x, 1.0); },
        max1_val, dmax1_dx, grid_sym
    );

    test_unary_F_compare(
        "min(x,1.0)",
        [](const FADF& x) { return fadbad::min(x, 1.0); },
        min1_val, dmin1_dx, grid_sym
    );
    
    test_unary_F_compare(
        "comp1(x) = erf(x)*xlog(x+2) + 0.5*abs(x-0.3)",
        [](const FADF& x) {
            return fadbad::erf(x) * fadbad::xlog(x + 2.0)
                 + 0.5 * fadbad::fabs(x - 0.3);
        },
        comp1_val, comp1_der, grid_comp_u1
    );

    test_unary_F_compare(
        "comp2(x) = max(x,1)*min(x+0.5,1) + x*x",
        [](const FADF& x) {
            return fadbad::max(x, 1.0) * fadbad::min(x + 0.5, 1.0)
                 + x * x;
        },
        comp2_val, comp2_der, grid_comp_u2
    );

    test_binary_F_compare(
        "comp3(x,y) = max(x,y)*(x+y) + min(x,1)*y^2 + erf(x-y)",
        [](const FADF& x, const FADF& y) {
            return fadbad::max(x, y) * (x + y)
                 + fadbad::min(x, 1.0) * y * y
                 + fadbad::erf(x - y);
        },
        comp3_val, comp3_dx, comp3_dy, grid_comp_x, grid_comp_y
    );

    test_binary_F_compare(
        "comp4(x,y) = abs(x-0.5)*max(y,1) + xlog(y+2.5) + x*y",
        [](const FADF& x, const FADF& y) {
            return fadbad::fabs(x - 0.5) * fadbad::max(y, 1.0)
                 + fadbad::xlog(y + 2.5)
                 + x * y;
        },
        comp4_val, comp4_dx, comp4_dy, grid_comp_x2, grid_comp_y2
    );
    
    // -----------------------------
    // B<double>
    // -----------------------------
    test_unary_B_compare(
        "fabs(x)",
        [](const BADB& x) { return fadbad::fabs(x); },
        abs_val, abs_der, grid_sym
    );

    test_unary_B_compare(
        "erf(x)",
        [](const BADB& x) { return fadbad::erf(x); },
        erf_val, erf_der, grid_sym
    );

    test_unary_B_compare(
        "xlog(x)  [assumed x*log(x)]",
        [](const BADB& x) { return fadbad::xlog(x); },
        xlog_val, xlog_der, grid_pos
    );

    test_binary_B_compare(
        "max(x,y)",
        [](const BADB& x, const BADB& y) { return fadbad::max(x, y); },
        max_val, dmax_dx, dmax_dy, grid2, grid2
    );

    test_binary_B_compare(
        "min(x,y)",
        [](const BADB& x, const BADB& y) { return fadbad::min(x, y); },
        min_val, dmin_dx, dmin_dy, grid2, grid2
    );

    test_unary_B_compare(
        "max(x,1.0)",
        [](const BADB& x) { return fadbad::max(x, 1.0); },
        max1_val, dmax1_dx, grid_sym
    );

    test_unary_B_compare(
        "min(x,1.0)",
        [](const BADB& x) { return fadbad::min(x, 1.0); },
        min1_val, dmin1_dx, grid_sym
    );

    test_unary_B_compare(
        "comp1(x) = erf(x)*xlog(x+2) + 0.5*abs(x-0.3)",
        [](const BADB& x) {
            return fadbad::erf(x) * fadbad::xlog(x + 2.0)
                 + 0.5 * fadbad::fabs(x - 0.3);
        },
        comp1_val, comp1_der, grid_comp_u1
    );

    test_unary_B_compare(
        "comp2(x) = max(x,1)*min(x+0.5,1) + x*x",
        [](const BADB& x) {
            return fadbad::max(x, 1.0) * fadbad::min(x + 0.5, 1.0)
                 + x * x;
        },
        comp2_val, comp2_der, grid_comp_u2
    );

    test_binary_B_compare(
        "comp3(x,y) = max(x,y)*(x+y) + min(x,1)*y^2 + erf(x-y)",
        [](const BADB& x, const BADB& y) {
            return fadbad::max(x, y) * (x + y)
                 + fadbad::min(x, 1.0) * y * y
                 + fadbad::erf(x - y);
        },
        comp3_val, comp3_dx, comp3_dy, grid_comp_x, grid_comp_y
    );

    test_binary_B_compare(
        "comp4(x,y) = abs(x-0.5)*max(y,1) + xlog(y+2.5) + x*y",
        [](const BADB& x, const BADB& y) {
            return fadbad::fabs(x - 0.5) * fadbad::max(y, 1.0)
                 + fadbad::xlog(y + 2.5)
                 + x * y;
        },
        comp4_val, comp4_dx, comp4_dy, grid_comp_x2, grid_comp_y2
    );

    return 0;
}
