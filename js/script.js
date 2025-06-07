

let Dimensions = 0;

function setODEs() {
    let dimInput = document.getElementById("dimensionField");

    Dimensions = parseInt(dimInput.options[dimInput.selectedIndex].value);

    for (let i = 0; i < Dimensions; i++) {
        let div = document.createElement("div");
        div.setAttribute("class", "div-parameter-member");

        let inputODEEl = document.createElement("input");
        inputODEEl.setAttribute("type", "text");
        inputODEEl.setAttribute("name", `ODE${i}`);
        inputODEEl.setAttribute("id", `ODE${i}`);
        inputODEEl.setAttribute("value", `D[${i}] = `);

        div.appendChild(inputODEEl);
        document.getElementById("div-ODEInput").appendChild(div);
    }

    let initial_solution_value_h5El = document.getElementById("h5-initialSolutionValueLabel");
    initial_solution_value_h5El.innerHTML = `X<sub>n</sub>(t<sub>0</sub>)`;

    for (let i = 0; i < Dimensions; i++) {
        let div = document.createElement("div");
        div.setAttribute("class", "div-parameter-member");

        let initial_value_inputEl = document.createElement("input");      
        initial_value_inputEl.setAttribute("type", "number");
        initial_value_inputEl.setAttribute("name", `Xn${i}`);
        initial_value_inputEl.setAttribute("id", `Xn${i}`);
        initial_value_inputEl.setAttribute("placeholder", `value of Xn[${i}] at t0`);

        div.appendChild(initial_value_inputEl);
        document.getElementById("div-solutionValue").appendChild(div);
    }

    document.getElementById("btnSetODEs").disabled = true;
    dimInput.disabled = true;
    document.getElementById("calcButton").disabled = false;
}

function calculate() {

    let eqn_input = "";

    for (let i = 0; i < Dimensions; i++) {
        let eqnInput = document.getElementById(`ODE${i}`);
        if (eqnInput.value.length == "") {
            alert(`Equation ${i+1} is missing.`);
            return;
        } else {
            eqn_input += eqnInput.value;
        }
    }

    let absTolInput = document.getElementById("abstol");
    let relTolInput = document.getElementById("reltol");
    let tspMethodInput = document.getElementById("timesteppingMethod");
    let t0ElVal = document.getElementById("init_t0value");
    let tfElVal = document.getElementById("init_tFinalvalue");

    let abstol = parseFloat(absTolInput.options[absTolInput.selectedIndex].value);
    let reltol = parseFloat(relTolInput.options[relTolInput.selectedIndex].value);
    let tspMethod = parseInt(tspMethodInput.options[tspMethodInput.selectedIndex].value);
    let t0 = parseFloat(t0ElVal.value);
    let tf = parseFloat(tfElVal.value);


    const info = {
        eqns: eqn_input,
        y0: new Array(Dimensions),
        T0: t0,
        Tf: tf, 
        absTol: abstol, 
        relTol: reltol, 
        timesteppingMethod: tspMethod
    };

    // Get Xn(t0)
    for (let i = 0; i < Dimensions; i++) {
        let X0 = document.getElementById(`Xn${i}`);
        info.y0[i] = parseFloat(X0.value);
    }

    // start our calculate animation
    let statusMessage = document.getElementById("statusMessage");
    statusMessage.innerText = "calculating .....";

    const workerSolver = new Worker('./js/workerSolver.js');

    workerSolver.postMessage(info);

    workerSolver.onmessage = function(e) {
        statusMessage.innerText = "";

        if (e.data.error !== undefined) {
            console.log(e.data.error);
            document.getElementById("errorReport").innerText = e.data.error;
        }

        if (e.data.t !== undefined && e.data.sol !== undefined) {
            if (!document.body.contains(document.getElementById("resultTable"))) {  // generate only one table
                genResultTable(e.data.t, e.data.sol);
            } else {
                document.getElementById("resultTable").remove();
                genResultTable(e.data.t, e.data.sol);
            }
        }
        workerSolver.terminate();
    };
}

function genResultTable(tData, solData) {
    
    // creates a <table> element and a <tbody> element
    const tbl = document.createElement("table");
    tbl.setAttribute("id", "resultTable");

    const tblBody = document.createElement("tbody");

    const headerText = ["Steps", "time (t)"];
    for (let i = 0; i < Dimensions; i++) {
        headerText.push(`Xn[${i}]`);
    }

    // create <th>
    const rowHeader = document.createElement("tr");

    for (let i = 0; i < headerText.length; i++) {
        const header = document.createElement("th");
        const ht = document.createTextNode(headerText[i]);
        header.appendChild(ht);
        rowHeader.appendChild(header);
    }

    tblBody.appendChild(rowHeader);

    // creating all cells
    for (let i = 0; i < solData.length; i++) {

        const row = document.createElement("tr");

        const stepCell = document.createElement("td");
        const stepCellText = document.createTextNode(`${i}`);
        stepCell.appendChild(stepCellText);
        row.appendChild(stepCell);

        const tCell = document.createElement("td");
        const tCellText = document.createTextNode(`${tData[i]}`);
        tCell.appendChild(tCellText);
        row.appendChild(tCell);

        for (let j = 0; j < Dimensions; j++) {
            // Create a <td> element and a text node, make the text
            // node the contents of the <td>, and put the <td> at
            // the end of the table row
            const cell = document.createElement("td");
            const cellText = document.createTextNode(`${solData[i][j]}`);
            cell.appendChild(cellText);
            row.appendChild(cell);
        }

        tblBody.appendChild(row);
    }

    // put the <tbody> in the <table>
    tbl.appendChild(tblBody);
    tbl.setAttribute("border", "2");
    tbl.setAttribute("align", "center");
    tbl.setAttribute("width", "80%");
    document.getElementById("div-resultTable").appendChild(tbl);
}
