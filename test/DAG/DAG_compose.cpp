#include <fstream>
#include <iomanip>

#include "ffexpr.hpp"
#include "ffunc.hpp"

///////////////////////////////////////////////////////////////////////////////

int
test_compose1()
{
  std::cout
      << "\n==============================================\ntest_compose1:\n";
  std::cout << "Compose G(X,Y) = sqr(Y)+exp(X) with Y <- F(X) = exp(X)\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");

  mc::FFVar F = exp(X);
  mc::FFVar G = sqr(Y) + F;

  auto sgF = DAG.subgraph({F});
  auto sgG = DAG.subgraph({G});
  DAG.output(sgG, " G");

  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  auto exG = mc::FFExpr::subgraph(&DAG, sgG);
  std::cout << "F(X) = " << exF[0] << "\nG(X,Y) = " << exG[0] << std::endl;

  // Compose G with Y <- F using array overload
  mc::FFVar const* pGoF = DAG.compose(1, &G, 1, &Y, &F);

  auto sgGoF = DAG.subgraph(1, pGoF);
  DAG.output(sgGoF, " GoF");
  auto exGoF = mc::FFExpr::subgraph(&DAG, sgGoF);
  std::cout << "GoF = " << exGoF[0] << "  (expected: sqr(exp(X))+exp(X))"
            << std::endl;

  delete[] pGoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_compose2()
{
  std::cout
      << "\n==============================================\ntest_compose2:\n";
  std::cout << "Compose G(Y,Z) = sqr(Y)+Z with Y <- F(X) and Z <- F(X)\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, Z;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");
  Z.set(&DAG, "Z");

  mc::FFVar F = exp(X);
  mc::FFVar G = sqr(Y) + Z;

  auto sgF = DAG.subgraph({F});
  auto sgG = DAG.subgraph({G});
  DAG.output(sgG, " G");

  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  auto exG = mc::FFExpr::subgraph(&DAG, sgG);
  std::cout << "F(X) = " << exF[0] << "\nG(Y,Z) = " << exG[0] << std::endl;

  // Compose using vector overload: Y <- F, Z <- F
  auto vGoF = DAG.compose({G}, {Y, Z}, {F, F});

  auto sgGoF = DAG.subgraph(vGoF);
  DAG.output(sgGoF, " GoF (vector)");
  auto exGoF = mc::FFExpr::subgraph(&DAG, sgGoF);
  std::cout << "GoF (vector) = " << exGoF[0] << std::endl;

  // Compose using array overload with variadic args: Y <- F, then Z <- F
  mc::FFVar const* pGoF = DAG.compose(1, &G, 1, &Y, &F, 1, &Z, &F);

  sgGoF = DAG.subgraph(1, pGoF);
  DAG.output(sgGoF, " GoF (array variadic)");
  exGoF = mc::FFExpr::subgraph(&DAG, sgGoF);
  std::cout << "GoF (array variadic) = " << exGoF[0]
            << "  (expected: sqr(exp(X))+exp(X))" << std::endl;

  delete[] pGoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_compose3()
{
  std::cout
      << "\n==============================================\ntest_compose3:\n";
  std::cout << "Compose with shared variables (original mov() bug case)\n";

  mc::FFGraph DAG;
  mc::FFVar x5, x30, x76, y5, y30, y76;
  x5.set(&DAG, "x5");
  x30.set(&DAG, "x30");
  x76.set(&DAG, "x76");
  y5.set(&DAG, "y5");
  y30.set(&DAG, "y30");
  y76.set(&DAG, "y76");

  mc::FFVar F = (-1.2040939193257074e+02) * (x5 * x76) +
                1.2040939193257074e+02 * (x5 * x30) +
                1.2040939193257074e+02 * x76;

  auto sgF = DAG.subgraph({F});
  DAG.output(sgF, " F");
  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  std::cout << "F = " << exF[0] << std::endl;

  // Compose: x5 <- y5, x30 <- y30, x76 <- y76
  auto F2 = DAG.compose({F}, {x5, x30, x76}, {y5, y30, y76});

  auto sgF2 = DAG.subgraph({F2});
  DAG.output(sgF2, " F2");
  auto exF2 = mc::FFExpr::subgraph(&DAG, sgF2);
  std::cout << "F2 = " << exF2[0]
            << "\n(F2 should be the same expression with y's replacing x's)"
            << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_insert()
{
  std::cout
      << "\n==============================================\ntest_insert:\n";
  std::cout << "Insert subgraphs from one DAG into another\n";

  mc::FFGraph DAG1;

  size_t const NX = 3;
  std::vector<mc::FFVar> X1(NX);
  for (auto& X1i : X1) X1i.set(&DAG1);

  std::vector<mc::FFVar> F1{exp(X1[1]), sqr(X1[0]) - exp(X1[1]) * X1[0],
                            sqr(X1[2]) - exp(X1[1]) * X1[2]};

  DAG1.output(DAG1.subgraph(F1), " F1");

  // Insert single dependents one at a time
  mc::FFGraph DAG2;
  std::vector<mc::FFVar> F2(2);
  DAG2.insert(&DAG1, 1, &F1[0], &F2[0]);
  DAG2.insert(&DAG1, 1, &F1[1], &F2[1]);
  DAG2.output(DAG2.subgraph(F2), " F2 (inserted one by one)");

  // Insert all dependents at once
  mc::FFGraph DAG3;
  std::vector<mc::FFVar> F3(3);
  DAG3.insert(&DAG1, F1, F3);
  DAG3.output(DAG3.subgraph(F3), " F3 (inserted all at once)");

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_substitute1()
{
  std::cout << "\n==============================================\ntest_"
               "substitute1:\n";
  std::cout << "Substitute a single AUX variable (intermediate result)\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, Z;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");
  Z.set(&DAG, "Z");

  // G = sqr(Y) + exp(X)  — exp(X) is an AUX node
  mc::FFVar eX = exp(X);
  mc::FFVar G  = sqr(Y) + eX;

  auto sgG = DAG.subgraph({G});
  DAG.output(sgG, " G");
  auto exG = mc::FFExpr::subgraph(&DAG, sgG);
  std::cout << "G = " << exG[0] << std::endl;

  // Substitute the AUX node exp(X) -> Z
  auto G2 = DAG.substitute({G}, {eX}, {Z});

  auto sgG2 = DAG.subgraph({G2});
  DAG.output(sgG2, " G2 = G with exp(X)->Z");
  auto exG2 = mc::FFExpr::subgraph(&DAG, sgG2);
  std::cout << "G2 = " << exG2[0] << "  (expected: sqr(Y)+Z)" << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_substitute2()
{
  std::cout << "\n==============================================\ntest_"
               "substitute2:\n";
  std::cout << "Substitute multiple AUX variables simultaneously\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, U, V;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");
  U.set(&DAG, "U");
  V.set(&DAG, "V");

  // Build F = sqr(X) * exp(Y) — both sqr(X) and exp(Y) are AUX nodes
  mc::FFVar sX = sqr(X);
  mc::FFVar eY = exp(Y);
  mc::FFVar F  = sX * eY;

  auto sgF = DAG.subgraph({F});
  DAG.output(sgF, " F");
  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  std::cout << "F = " << exF[0] << std::endl;

  // Substitute both AUX nodes: sqr(X) -> U, exp(Y) -> V
  auto F2 = DAG.substitute({F}, {sX, eY}, {U, V});

  auto sgF2 = DAG.subgraph({F2});
  DAG.output(sgF2, " F2 = F with sqr(X)->U, exp(Y)->V");
  auto exF2 = mc::FFExpr::subgraph(&DAG, sgF2);
  std::cout << "F2 = " << exF2[0] << "  (expected: U*V)" << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_substitute3()
{
  std::cout << "\n==============================================\ntest_"
               "substitute3:\n";
  std::cout
      << "Substitute VAR variables (substitute works for both VAR and AUX)\n";

  // Repeat the compose4 test case using substitute instead of compose
  mc::FFGraph DAG;
  mc::FFVar x5, x30, x76, y5, y30, y76;
  x5.set(&DAG, "x5");
  x30.set(&DAG, "x30");
  x76.set(&DAG, "x76");
  y5.set(&DAG, "y5");
  y30.set(&DAG, "y30");
  y76.set(&DAG, "y76");

  mc::FFVar F = (-1.2040939193257074e+02) * (x5 * x76) +
                1.2040939193257074e+02 * (x5 * x30) +
                1.2040939193257074e+02 * x76;

  auto sgF = DAG.subgraph({F});
  DAG.output(sgF, " F");
  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  std::cout << "F = " << exF[0] << std::endl;

  auto F2 = DAG.substitute({F}, {x5, x30, x76}, {y5, y30, y76});

  auto sgF2 = DAG.subgraph({F2});
  DAG.output(sgF2, " F2 = F with x5->y5, x30->y30, x76->y76");
  auto exF2 = mc::FFExpr::subgraph(&DAG, sgF2);
  std::cout << "F2 = " << exF2[0]
            << "\n(F2 should be the same expression with y's replacing x's)"
            << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_substitute4()
{
  std::cout << "\n==============================================\ntest_"
               "substitute4:\n";
  std::cout << "Substitute an AUX node used multiple times (pruned subgraph)\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, Z;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");
  Z.set(&DAG, "Z");

  // t = X*Y is used twice: F = t + sqr(t)
  mc::FFVar t = X * Y;
  mc::FFVar F = t + sqr(t);

  auto sgF = DAG.subgraph({F});
  DAG.output(sgF, " F");
  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  std::cout << "F = " << exF[0] << std::endl;

  // Substitute t -> Z  (the AUX node X*Y is replaced by the VAR Z)
  auto F2 = DAG.substitute({F}, {t}, {Z});

  auto sgF2 = DAG.subgraph({F2});
  DAG.output(sgF2, " F2 = F with X*Y->Z");
  auto exF2 = mc::FFExpr::subgraph(&DAG, sgF2);
  std::cout << "F2 = " << exF2[0] << "  (expected: Z+sqr(Z))" << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_substitute5()
{
  std::cout << "\n==============================================\ntest_"
               "substitute5:\n";
  std::cout << "Substitute with variadic args (multiple substitution groups)\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, U, V;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");
  U.set(&DAG, "U");
  V.set(&DAG, "V");

  // F = exp(X) + log(Y) — both are AUX nodes
  mc::FFVar eX = exp(X);
  mc::FFVar lY = log(Y);
  mc::FFVar F  = eX + lY;

  auto sgF = DAG.subgraph({F});
  DAG.output(sgF, " F");
  auto exF = mc::FFExpr::subgraph(&DAG, sgF);
  std::cout << "F = " << exF[0] << std::endl;

  // Substitute using two separate groups via variadic args:
  //   group 1: exp(X) -> U
  //   group 2: log(Y) -> V
  // Note: variadic args require explicit std::vector (initializer lists
  // cannot be deduced as template arguments in a parameter pack)
  std::vector<mc::FFVar> vAuxTarg2{lY}, vAuxSubst2{V};
  auto F2 = DAG.substitute({F}, {eX}, {U}, vAuxTarg2, vAuxSubst2);

  auto sgF2 = DAG.subgraph(F2);
  DAG.output(sgF2, " F2 = F with exp(X)->U, log(Y)->V (variadic)");
  auto exF2 = mc::FFExpr::subgraph(&DAG, sgF2);
  std::cout << "F2 = " << exF2[0] << "  (expected: U+V)" << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
test_substitute6()
{
  std::cout << "\n==============================================\ntest_"
               "substitute6:\n";
  std::cout << "Verify compose rejects AUX targets with NOTVAR exception\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, Z;
  X.set(&DAG, "X");
  Y.set(&DAG, "Y");
  Z.set(&DAG, "Z");

  mc::FFVar eX = exp(X);
  mc::FFVar G  = sqr(Y) + eX;

  // Attempt to compose with an AUX target — should throw NOTVAR
  try
  {
    auto G2 = DAG.compose({G}, {eX}, {Z});
    std::cout << "ERROR: compose did NOT throw for AUX target!" << std::endl;
    return 1;
  }
  catch (mc::FFBase::Exceptions& eObj)
  {
    if (eObj.ierr() == mc::FFBase::Exceptions::NOTVAR)
    {
      std::cout << "PASS: compose correctly threw NOTVAR: " << eObj.what()
                << std::endl;
    }
    else
    {
      std::cout << "FAIL: compose threw unexpected exception: " << eObj.what()
                << std::endl;
      return 1;
    }
  }

  // Same substitution via substitute() should succeed
  auto G2   = DAG.substitute({G}, {eX}, {Z});
  auto sgG2 = DAG.subgraph({G2});
  DAG.output(sgG2, " G2 via substitute");
  auto exG2 = mc::FFExpr::subgraph(&DAG, sgG2);
  std::cout << "PASS: substitute succeeded: G2 = " << exG2[0] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main()
{
  bool failed = true;

  try
  {
    test_compose1();
    test_compose2();
    test_compose3();
    test_insert();
    test_substitute1();
    test_substitute2();
    test_substitute3();
    test_substitute4();
    test_substitute5();
    test_substitute6();
    failed = false;
  }
  catch (mc::FFBase::Exceptions& eObj)
  {
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborting." << std::endl;
  }
  catch (...)
  {
    std::cerr << "Error during DAG evaluation\n"
              << "Aborting." << std::endl;
  }

  std::cout << "\n=== Results: " << (failed ? "failed" : "passed") << " ===\n";
  return failed;
}
