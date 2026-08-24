/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "MC++", "index.html", [
    [ "MC++ (version 3.0): Toolkit for Construction, Manipulation and Bounding of Factorable Functions", "index.html", "index" ],
    [ "Non-Verified Interval Arithmetic for Factorable Functions", "page_INTERVAL.html", [
      [ "How do I compute interval bounds on the range of a factorable function?", "page_INTERVAL.html#sec_INTERVAL_use", null ],
      [ "Which functions are overloaded in mc::Interval?", "page_INTERVAL.html#sec_INTERVAL_fct", null ],
      [ "What are the options in mc::Interval and how are they set?", "page_INTERVAL.html#sec_INTERVAL_opt", null ],
      [ "What errors can be encountered in using mc::Interval?", "page_INTERVAL.html#sec_INTERVAL_err", null ],
      [ "References", "page_INTERVAL.html#sec_INTERVAL_refs", null ]
    ] ],
    [ "McCormick Relaxation Arithmetic for Factorable Functions", "page_MCCORMICK.html", [
      [ "How do I compute McCormick relaxations of a factorable function?", "page_MCCORMICK.html#sec_MCCORMICK_use", null ],
      [ "How do I compute a subgradient of the McCormick relaxations?", "page_MCCORMICK.html#sec_MCCORMICK_sub", null ],
      [ "How do I compute McCormick relaxations of the partial derivatives or the Taylor coefficients of a factorable function using FADBAD++?", "page_MCCORMICK.html#sec_MCCORMICK_fadbad", null ],
      [ "Which functions are overloaded in McCormick relaxation arithmetic?", "page_MCCORMICK.html#sec_MCCORMICK_fct", null ],
      [ "What are the options in mc::McCormick and how are they set?", "page_MCCORMICK.html#sec_MCCORMICK_opt", null ],
      [ "What Errors Can Be Encountered during the Computation of Convex/Concave Bounds?", "page_MCCORMICK.html#sec_MC_err", null ],
      [ "References", "page_MCCORMICK.html#sec_MC_refs", null ]
    ] ],
    [ "Eigenvalue Arithmetic for Factorable Functions", "page_SPECBND.html", [
      [ "How do I compute (interval) bounds on the spectrum of the Hessian matrix of a factorable function?", "page_SPECBND.html#sec_SPECBND_I", null ],
      [ "How do I compute convex/concave relaxations on the spectrum of the Hessian matrix of a factorable function?", "page_SPECBND.html#sec_SPECBND_MC", null ],
      [ "Which functions are overloaded in mc::Specbnd eigenvalue arithmetic?", "page_SPECBND.html#sec_SPECBND_fct", null ],
      [ "How do I compute spectral bounds of the partial derivatives or the Taylor coefficients of a factorable function using FADBAD++?", "page_SPECBND.html#sec_SPECBND_fadbad", null ],
      [ "How are the options set for the computation of a spectral bound?", "page_SPECBND.html#sec_SPECBND_opt", null ],
      [ "Errors What errors can I encounter during computation of a spectral bound?", "page_SPECBND.html#sec_SPECBND_err", null ],
      [ "References", "page_SPECBND.html#sec_SPECBND_refs", null ]
    ] ],
    [ "Taylor Model Arithmetic for Factorable Functions", "page_TAYLOR.html", [
      [ "How do I compute a Taylor model with interval remainder bound of a factorable function?", "page_TAYLOR.html#sec_TAYLOR_I", null ],
      [ "How do I compute a Taylor model with convex/concave remainder bounds of a factorable function?", "page_TAYLOR.html#sec_TAYLOR_MC", null ],
      [ "How do I compute Taylor models of the partial derivatives or the Taylor coefficients of a factorable function using FADBAD++?", "page_TAYLOR.html#sec_TAYLOR_fadbad", null ],
      [ "Which functions are overloaded for Taylor model arithmetic?", "page_TAYLOR.html#sec_TAYLOR_fct", null ],
      [ "How are the options set for the computation of a Taylor model?", "page_TAYLOR.html#sec_TAYLOR_opt", null ],
      [ "Errors What errors can I encounter during computation of a Taylor model?", "page_TAYLOR.html#sec_TM_err", null ],
      [ "References", "page_TAYLOR.html#sec_TM_refs", null ]
    ] ],
    [ "Chebyshev Model Arithmetic for Factorable Functions", "page_CHEBYSHEV.html", "page_CHEBYSHEV" ],
    [ "Ellipsoidal Calculus and Ellipsoidal Arithmetic for Factorable Functions", "page_ELLIPSOID.html", [
      [ "How do I define an ellipsoid and apply ellipsoidal calculus?", "page_ELLIPSOID.html#sec_ELLCALC", null ],
      [ "How do I compute an ellipsoidal enclosure for the image set of an ellipsoid under a factorable function?", "page_ELLIPSOID.html#sec_ELLIMG", null ],
      [ "How are the options set for ellipsoidal calculus and ellipsoidal arithmetic?", "page_ELLIPSOID.html#sec_ELL_opt", null ],
      [ "What are the errors encountered during ellipsoidal calculus and ellipsoidal arithmetic?", "page_ELLIPSOID.html#sec_ELL_err", null ],
      [ "References", "page_ELLIPSOID.html#sec_ELL_refs", null ]
    ] ],
    [ "Polyhedral Arithmetic for Factorable Functions", "page_POLYHEDRAL.html", [
      [ "What is the theory behind polyhedral relaxations?", "page_POLYHEDRAL.html#sec_POLIMG_THEOR", null ],
      [ "How to generate a polyhedral relaxation of a factorable function?", "page_POLYHEDRAL.html#sec_POLIMG_COMP", null ],
      [ "What are the options in computing a polyhedral relaxation?", "page_POLYHEDRAL.html#sec_POLIMG_OPT", null ],
      [ "What errors may be encoutered in computing a polyhedral relaxation?", "page_POLYHEDRAL.html#sec_POLIMG_ERR", null ],
      [ "References", "page_POLYHEDRAL.html#sec_POLIMG_REFS", null ]
    ] ],
    [ "Interval Superposition Model Arithmetic for Factorable Functions", "page_ISM.html", [
      [ "How do I compute an ISM of a factorable function?", "page_ISM.html#sec_ISM_use", null ],
      [ "Errors Errors encountered during computation of an ISM", "page_ISM.html#sec_ISM_err", null ],
      [ "References", "page_ISM.html#sec_ISM_refs", null ]
    ] ],
    [ "Dependence Structure Detection for Factorable Functions", "page_FFDEP.html", [
      [ "How Do I Determine the Structure of a Factorable Function?", "page_FFDEP.html#sec_FFDepEval", null ],
      [ "Errors Encountered in Determining the Structure of a Factorable Function?", "page_FFDEP.html#sec_FFDepErr", null ]
    ] ],
    [ "Invertible Structure Detection in Factorable Functions", "page_FFINV.html", [
      [ "How Do I Determine the Invertibility of a Factorable Function?", "page_FFINV.html#sec_FFInvEval", null ],
      [ "Available Options for Invertible Structure Detection", "page_FFINV.html#sec_FFInvOpt", null ],
      [ "Errors Encountered during Invertible Structure Detection", "page_FFINV.html#sec_FFInvErr", null ]
    ] ],
    [ "Construction, Manipulation and Evaluation of Factorable Functions", "page_FFUNC.html", "page_FFUNC" ],
    [ "Bug List", "bug.html", null ],
    [ "Topics", "topics.html", "topics" ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Enumerator", "functions_eval.html", "functions_eval" ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"classmc_1_1FFGraph.html#a9103357bca78683c43e72e644474ef67",
"classmc_1_1PolCut.html#aec222333ec4b858417174aeb08d89ba4a5bbe922e44234866dd9ab192681d1943",
"classmc_1_1SCModel.html#ab583baa0a29f6d382f389d737671319f",
"classmc_1_1SPoly.html#a666104111215753bfcca03cee4fac686",
"functions_func_b.html",
"group__FFunc.html#ga84f601189c93e47bb32f074f341a1610",
"group__SICHEBYSHEV.html#gae0e0ce3fdfbbfa6d19721db08508c2e9",
"structmc_1_1Ellipsoid_1_1Options.html#abcaa7e1f1672cb03e3aea0813e745297",
"structmc_1_1SQuad_1_1Options.html#af223be48b0841a8a163bf9de63398fd2"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';