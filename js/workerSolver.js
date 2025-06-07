
let Dimensions;

const Dopri54 = {
    c2: 1.0 / 5.0,
    c3: 3.0 / 10.0,
    c4: 4.0 / 5.0,
    c5: 8.0 / 9.0,
    a21: 1.0 / 5.0,
    a31: 3.0 / 40.0,
    a32: 9.0 / 40.0,
    a41: 44.0 / 45.0,
    a42: -56.0 / 15.0,
    a43: 32.0 / 9.0,
    a51: 19372.0 / 6561.0,
    a52: -25360.0 / 2187.0,
    a53: 64448.0 / 6561.0,
    a54: -212.0 / 729.0,
    a61: 9017.0 / 3168.0,
    a62: -355.0 / 33.0,
    a63: 46732.0 / 5247.0,
    a64: 49.0 / 176.0,
    a65: -5103.0 / 18656.0,
    a71: 35.0 / 384.0,
    a73: 500.0 / 1113.0,
    a74: 125.0 / 192.0,
    a75: -2187.0 / 6784.0,
    a76: 11.0 / 84.0,
    b1: 35.0 / 384.0,
    b3: 500.0 / 1113.0,
    b4: 125.0 / 192.0,
    b5: -2187.0 / 6784.0,
    b6: 11.0 / 84.0,
    e1: 5179.0 / 57600.0,
    e3: 7571.0 / 16695.0,
    e4: 393.0 / 640.0,
    e5: -92097.0 / 339200.0,
    e6: 187.0 / 2100.0,
    e7: 1.0 / 40.0
};


const Controllers = Object.freeze({
    STANDARD: 0,
    H211PI: 1,
    H312PID: 2,
    H211B: 3,
    H312B: 4,
    PI42: 5,
    H0321: 6,
    H321: 7


    /*
    ***Recommended Controllers with Stepsize Low-Pass Filters and their Problem Classes***
    *  Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003)
    *  https://dl.acm.org/doi/10.1145/641876.641877 -> Table III. page 24
    */

    /*--------------------------------------------------------------------------
    * kbeta1 | kbeta2 | kbeta3 | alpha2 | alpha3 | Class    | Problem Type
    *-------------------------------------------------------------------------
    * 1/b    | 1/b    | 0      | 1/b    | 0      | H211b    | medium to nonsmooth
    * 1/6    | 1/6    | 0      | 0      | 0      | H211 PI  | medium to nonsmooth
    * 1/b    | 2/b    | 1/b    | 3/b    | 1/b    | H312b    | nonsmooth
    * 1/18   | 1/9    | 1/18   | 0      | 0      | H312 PID | nonsmooth
    * 5/4    | 1/2    |-3/4    |-1/4    |-3/4    | H0321    | smooth
    * 1/3    | 1/18   |-5/18   |-5/6    |-1/6    | H321     | medium
    *-------------------------------------------------------------------------
    */
});


class Solver {
    #m_p = 4.0;                  // the order corresponding to the RK method
    #m_k = this.#m_p + 1.0;      // EPS => k = p + 1 and EPUS => k = p
    #m_kappa = 1.5;              // kappa ∈ [0.7, 2] as suggested in the literature
    #m_acceptSF = 0.81;          // accept safety factor

    #m_h;
    #m_t;
    #m_tFinal;
    #m_absTol;
    #m_relTol;
    #m_beta1;
    #m_beta2;
    #m_beta3;
    #m_alpha2;
    #m_alpha3;
    #m_denseOut;

    #m_cerr1 = 1.0;
    #m_cerr2 = 1.0;
    #m_rh1 = 1.0;
    #m_rh2 = 1.0;

    #m_rejectedSteps = 0;
    #m_acceptedSteps = 0;

    #m_yn = new Array(Dimensions);
    #m_X = new Array(Dimensions);
    #m_K1 = new Array(Dimensions);
    #m_K2 = new Array(Dimensions);
    #m_K3 = new Array(Dimensions);
    #m_K4 = new Array(Dimensions);
    #m_K5 = new Array(Dimensions);
    #m_K6 = new Array(Dimensions);
    #m_K7 = new Array(Dimensions);
    #m_ynew = new Array(Dimensions);
    #m_sci = new Array(Dimensions);
    #m_truncationErrors = new Array(Dimensions);

    #inputEquationString;

    #terminateSolver = false;


    #m_F(t, Xn, D) {

        try {
            eval(this.#inputEquationString);
        } catch (error) {
            this.#terminateSolver = true;
            postMessage({ error: error.message });
        }
    }

    #set_controller_parameters(b1, b2, b3, a2, a3) {
        this.#m_beta1 = b1 / this.#m_k;
        this.#m_beta2 = b2 / this.#m_k;
        this.#m_beta3 = b3 / this.#m_k;
        this.#m_alpha2 = a2;
        this.#m_alpha3 = a3;
    }

    #hairer_norm(a, sci) {
        let sum_of_sqrd = 0.0;
        for (let i = 0; i < Dimensions; i++) {
            sum_of_sqrd = sum_of_sqrd + Math.pow(a[i] / sci[i], 2.0);
        }
        return Math.sqrt(sum_of_sqrd / Dimensions);
    }

    #initialize_stepsize(t0, y0)
    {
        /*
         * (a) Do one function evaluation f(t0, y0) at the initial point.
         * Then put d0 = ||y0|| and d1 = ||f(t0, y0)||, using the hairer's norm
         * with sci = Atol + |y0_i| * Rtol
        */
        let f1 = new Array(Dimensions);
        this.#m_F(t0, y0, f1);

        let sci = new Array(Dimensions);
        for (let i = 0; i < Dimensions; i++)
        {
            sci[i] = this.#m_absTol + Math.abs(y0[i]) * this.#m_relTol;
        }

        let d0 = this.#hairer_norm(y0, sci);
        let d1 = this.#hairer_norm(f1, sci);

        // (b) As a first guess for the step size let
        let h0 = 0.01 * (d0 / d1);

        // If either d0 or d1 is < 1e-5 we put h0 = 1e-6
        if (d0 < 1e-5 || d1 < 1e-5)
        {
            h0 = 1e-6;
        }

        // (c) Perform one explicit Euler stpe, y1 = y0 + h0 * f(t0, y0)
        let y1 = new Array(Dimensions);
        for (let i = 0; i < Dimensions; i++)
        {
            y1[i] = y0[i] + h0 * f1[i];
        }

        let f2 = new Array(Dimensions);;
        this.#m_F(t0 + h0, y1, f2);

        /*
         * (d) Compute d2 = ||f(t0 + h0, y1) - f(t0, y0)|| / h0 as an
         * estimate of the second derivative of the solution; 
         * again by using hairer's norm.
        */
        let diff_f2f1 = new Array(Dimensions);;
        for (let i = 0; i < Dimensions; i++)
        {
            diff_f2f1[i] = Math.abs(f2[i] - f1[i]);
        }

        let d2 = this.#hairer_norm(diff_f2f1, sci) / h0;

        let max_d1d2 = Math.max(d1, d2);

        /*
         * (e) Compute a step size h1 from the relation, h1^(p+1) * max(d1, d2) = 0.01, 
         * where p - is order of the method. If max(d1, d2) <= 10^-15, 
         * put h1 = max(10^-6, h0 * 10^-3);
        */
        let h1 = Math.pow(10.0, (-2.0 - Math.log10(max_d1d2)) / this.#m_k);
        if (max_d1d2 <= 1e-15)
        {
            h1 = Math.max(1e-6, h0 * 1e-3);
        }

        // f. Finally we propose as starting step size
        let h = Math.min(100.0 * h0, h1);

        return h;
    }

    #process() {

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i];
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t, this.#m_X, this.#m_K1); //--------------------- 1ST-stage -----------------------
        }

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i] + this.#m_h * (Dopri54.a21 * this.#m_K1[i]);
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t + Dopri54.c2 * this.#m_h, this.#m_X, this.#m_K2); //--------------------- 2ND-stage -----------------------
        }

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i] + this.#m_h * (Dopri54.a31 * this.#m_K1[i] + Dopri54.a32 * this.#m_K2[i]);
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t + Dopri54.c3 * this.#m_h, this.#m_X, this.#m_K3); //--------------------- 3RD-stage -----------------------
        }

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i] + this.#m_h * (Dopri54.a41 * this.#m_K1[i] + Dopri54.a42 * this.#m_K2[i] +
                Dopri54.a43 * this.#m_K3[i]);
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t + Dopri54.c4 * this.#m_h, this.#m_X, this.#m_K4); //--------------------- 4TH-stage -----------------------
        }

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i] + this.#m_h * (Dopri54.a51 * this.#m_K1[i] + Dopri54.a52 * this.#m_K2[i] +
                Dopri54.a53 * this.#m_K3[i] + Dopri54.a54 * this.#m_K4[i]);
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t + Dopri54.c5 * this.#m_h, this.#m_X, this.#m_K5); //--------------------- 5TH-stage -----------------------
        }

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i] + this.#m_h * (Dopri54.a61 * this.#m_K1[i] + Dopri54.a62 * this.#m_K2[i] +
                Dopri54.a63 * this.#m_K3[i] + Dopri54.a64 * this.#m_K4[i] + Dopri54.a65 * this.#m_K5[i]);
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t + this.#m_h, this.#m_X, this.#m_K6); //--------------------- 6TH-stage -----------------------
        }

        for (let i = 0; i < Dimensions; i++) {
            this.#m_X[i] = this.#m_yn[i] + this.#m_h * (Dopri54.a71 * this.#m_K1[i] + Dopri54.a73 * this.#m_K3[i] +
                Dopri54.a74 * this.#m_K4[i] + Dopri54.a75 * this.#m_K5[i] + Dopri54.a76 * this.#m_K6[i]);
        }

        // Calculate the 5th-order and 4th-order accurate solution.
        for (let i = 0; i < Dimensions; i++) {
            this.#m_ynew[i] = this.#m_yn[i] + this.#m_h * (Dopri54.b1 * this.#m_K1[i] + Dopri54.b3 * this.#m_K3[i] +
                Dopri54.b4 * this.#m_K4[i] + Dopri54.b5 * this.#m_K5[i] + Dopri54.b6 * this.#m_K6[i]); // 5th-Order accurate solution. Used to advance the solution.
        }

        if (!this.#terminateSolver) {
            this.#m_F(this.#m_t + this.#m_h, this.#m_ynew, this.#m_K7); //--------------------- 7TH-stage ----------------------- reuse for error-estimation
        }


        // Calculate local errors
        for (let i = 0; i < Dimensions; i++) {
            this.#m_truncationErrors[i] = this.#m_h * ((Dopri54.b1 - Dopri54.e1) * this.#m_K1[i] + (Dopri54.b3 - Dopri54.e3) * this.#m_K3[i] +
                                                       (Dopri54.b4 - Dopri54.e4) * this.#m_K4[i] + (Dopri54.b5 - Dopri54.e5) * this.#m_K5[i] +
                                                       (Dopri54.b6 - Dopri54.e6) * this.#m_K6[i] - Dopri54.e7 * this.#m_K7[i]);
        }

        for (let i = 0; i < Dimensions; i++) {
            // absTol and relTol are the desired tolerances prescribed by the user.
            this.#m_sci[i] = this.#m_absTol + Math.max(Math.abs(this.#m_yn[i]), Math.abs(this.#m_ynew[i])) * this.#m_relTol;
        }

        // local error is controlled by error-per-unit-steps (EPUS) or error-per-steps (EPS)
        // let ees = this.#hairer_norm() / this.#m_h;  // Error-per-unit-steps (EPUS)
        let ees = this.#hairer_norm(this.#m_truncationErrors, this.#m_sci); // Error-per-steps (EPS)

        return ees;
    }

    #filter(cerrPres, cerrOld1, cerrOld2, rho1, rho2) {
        /**
         * The General controller formula for Order-Dynamics pD <= 3 with Control error filtering.
         * Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003)
         * https://dl.acm.org/doi/10.1145/641876.641877 -> page 22
         */
        let result = Math.pow(cerrPres, this.#m_beta1) *
            Math.pow(cerrOld1, this.#m_beta2) *
            Math.pow(cerrOld2, this.#m_beta3) *
            Math.pow(rho1, -this.#m_alpha2) *
            Math.pow(rho2, -this.#m_alpha3);
        return result;
    }

    #controlStepSize(ratio) {
        this.#m_h *= ratio;
    }

    #interpolate(theta, hPresent) {
        const C1 = 5.0 * (2558722523.0 - 31403016.0 * theta) / 11282082432.0;
        const C3 = 100.0 * (882725551.0 - 15701508.0 * theta) / 32700410799.0;
        const C4 = 25.0 * (443332067.0 - 31403016.0 * theta) / 1880347072.0;
        const C5 = 32805.0 * (23143187.0 - 3489224.0 * theta) / 199316789632.0;
        const C6 = 55.0 * (29972135.0 - 7076736.0 * theta) / 822651844.0;
        const C7 = 10.0 * (7414447.0 - 829305.0 * theta) / 29380423.0;

        let theta_sqr = Math.pow(theta, 2.0);
        let term1 = theta_sqr * (3.0 - 2.0 * theta);
        let term2 = theta_sqr * Math.pow(theta - 1.0, 2.0);
        let term3 = theta * Math.pow(theta - 1.0, 2.0);
        let term4 = (theta - 1.0) * Math.pow(theta, 2.0);

        let b1Theta = term1 * Dopri54.b1 + term3 - term2 * C1;
        let b3Theta = term1 * Dopri54.b3 + term2 * C3;
        let b4Theta = term1 * Dopri54.b4 - term2 * C4;
        let b5Theta = term1 * Dopri54.b5 + term2 * C5;
        let b6Theta = term1 * Dopri54.b6 - term2 * C6;
        let b7Theta = term4 + term2 * C7;

        let solution = new Array(Dimensions);
        for (let i = 0; i < Dimensions; i++) {
            /* code */
            solution[i] = this.#m_yn[i] + hPresent * (b1Theta * this.#m_K1[i] + b3Theta * this.#m_K3[i] + b4Theta * this.#m_K4[i] +
                b5Theta * this.#m_K5[i] + b6Theta * this.#m_K6[i] + b7Theta * this.#m_K7[i]);
        }

        return solution;
    }

    m_tOut = new Array();
    m_yOut = new Array();

    // solver(controllerType, f(t, y), y0, h, t0, tFinal, absTolerance=1E-6, relTolerance=1E-4, denseOut=false)
    constructor(controller, y0, t0, tFinal, absTol = 1e-6, relTol = 1e-4, denseOut = false, input_eqn_string) {
        this.#m_t = t0;
        this.#m_yn = y0;
        this.#m_tFinal = tFinal;
        this.#m_absTol = absTol;
        this.#m_relTol = relTol;
        this.#m_denseOut = denseOut;
        this.#inputEquationString = input_eqn_string;
        this.#m_h = this.#initialize_stepsize(t0, y0);

        switch (controller) {
            case 0:  // Standard
                this.#set_controller_parameters(1.0, 0.0, 0.0, 0.0, 0.0);
                break;
            case 1:  // H211PI
                this.#set_controller_parameters(1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0);
                break;
            case 2:  // H312PID
                this.#set_controller_parameters(1.0 / 18.0, 1.0 / 9.0, 1.0 / 18.0, 0.0, 0.0);
                break;
            case 3:  // H211B
                this.#set_controller_parameters(1.0 / 4.0, 1.0 / 4.0, 0.0, 1.0 / 4.0, 0.0);
                break;
            case 4:  // H312B
                this.#set_controller_parameters(1.0 / 8.0, 2.0 / 8.0, 1.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0);
                break;
            case 5:  // PI42
                this.#set_controller_parameters(3.0 / 5.0, -1.0 / 5.0, 0.0, 0.0, 0.0);
                break;
            case 6:  // H0321
                this.#set_controller_parameters(5.0 / 4.0, 1.0 / 2.0, -3.0 / 4.0, -1.0 / 4.0, -3.0 / 4.0);
                break;
            case 7:  // H321
                this.#set_controller_parameters(1.0 / 3.0, 1.0 / 18.0, -5.0 / 18.0, -5.0 / 6.0, -1.0 / 6.0);
            default:
                break;
        }
    }

    solve() {

        this.m_tOut.push(this.#m_t);
        let cpy = [...this.#m_yn];
        this.m_yOut.push(cpy);

        // Closed-loop system
        while (this.#m_t < this.#m_tFinal) {

            this.#m_h = Math.min(this.#m_h, this.#m_tFinal - this.#m_t);

            let err_estimate_scaled = this.#process();

            if (this.#terminateSolver) {
                break;
            }

            let cerr = 1.0 / err_estimate_scaled;
            let rho = this.#filter(cerr, this.#m_cerr1, this.#m_cerr2, this.#m_rh1, this.#m_rh2);

            // save previous values for the next step.
            this.#m_cerr2 = this.#m_cerr1;
            this.#m_cerr1 = cerr;

            this.#m_rh2 = this.#m_rh1;
            this.#m_rh1 = rho;

            // apply a limiter
            let ratio = 1.0 + this.#m_kappa * Math.atan((rho - 1.0) / this.#m_kappa);

            if (ratio < this.#m_acceptSF) {
                // reject steps and recalcualte with the new stepsize
                this.#controlStepSize(ratio);
                this.#m_rejectedSteps++;
                continue;
            } else {
                // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.
                if (this.#m_denseOut) {
                    // do 1-more extra steps at the middle between yn and yn1.
                    /*
                     * theta ∈ [0, 1], theta = 0 => yn, theta = 1 => yn1
                     * Interpolate at open-interval theta ∈ (0, 1)
                     * un+1(t + theta*h) = yn + h * sum(bi(theta)*Ki), i = 1...s, theta ∈ (0, 1)
                     */
                    let theta = 0.5;
                    let extraSteps = this.#interpolate(theta, this.#m_h);
                    this.m_tOut.push(this.#m_t + theta * this.#m_h);
                    this.m_yOut.push(extraSteps);
                }

                this.#m_t += this.#m_h;
                this.#controlStepSize(ratio);
                this.#m_acceptedSteps++;

                for (let i = 0; i < Dimensions; i++) {
                    this.#m_yn[i] = this.#m_ynew[i];
                }

                this.m_tOut.push(this.#m_t);
                let cpy = [...this.#m_ynew];
                this.m_yOut.push(cpy);
            }
        }
    }

    // displaySteps() {
    //     console.log(`\n\tsteps: accepted ${this.#m_acceptedSteps} rejected ${this.#m_rejectedSteps}`);
    // }

    // displayResults() {

    //     console.log('\n');

    //     for (let i = 0; i < this.m_yOut.length; i++) {
    //         console.log(`step ${i} at t = ${this.m_tOut[i]}\n`);
    //         for (let j = 0; j < Dimensions; j++) {
    //             console.log(`${this.m_yOut[i][j]}`);
    //         }
    //         console.log('\n');
    //     }
    // }
}

function startSolver(equation_input, y0, t0, tf, abstol, reltol, method) {

    Dimensions = y0.length;
    const solver = new Solver(method, y0, t0, tf, abstol, reltol, false, equation_input);

    solver.solve();

    return { tOut: solver.m_tOut, solOut: solver.m_yOut };
}

onmessage = function (e) {
    let equation_input = e.data.eqns;
    let y0 = e.data.y0;
    let t0 = e.data.T0;
    let tf = e.data.Tf;
    let abstol = e.data.absTol;
    let reltol = e.data.relTol;
    let method = e.data.timesteppingMethod;

    console.log(method);
    const solverRes = startSolver(equation_input, y0, t0, tf, abstol, reltol, method);
    postMessage({ t: solverRes.tOut, sol: solverRes.solOut });
};