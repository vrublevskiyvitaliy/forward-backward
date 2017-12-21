document.addEventListener('DOMContentLoaded',  () => {
    const initForm = document.getElementById('initForm');
    const matrixForm = document.getElementById('matrixForm');
    const resultBlock = document.getElementById('resultBlock');

    const fieldN = document.getElementById('varN'); // number of observations
    const fieldK = document.getElementById('varK'); // number of states
    const obsField = document.getElementById('obsField');

    const vectorPi = document.getElementById('vectorPi');
    const matrixA = document.getElementById('matrixA');
    const matrixB = document.getElementById('matrixB');

    let ObservationN, StatesN;

    initForm.addEventListener('submit', (event) => {
        event.preventDefault();

        const N = parseInt(fieldN.value);
        const K = parseInt(fieldK.value);

        ObservationN = parseInt(fieldN.value);
        StatesN = parseInt(fieldK.value);

        vectorPi.innerHTML = generateInputs(K);
        matrixA.innerHTML = generateMatrix(K, K);
        matrixB.innerHTML = generateMatrix(K, N);

        matrixForm.classList.add('is-visible');
        resultBlock.classList.remove('is-visible');
    });

    matrixForm.addEventListener('submit', (event) => {
        event.preventDefault();

        const result = Forward_Backward(
            getObservations(obsField),
            getInputValues(vectorPi),
            getMatrixValues(matrixA),
            getMatrixValues(matrixB)
        );

        resultBlock.innerText = result;
        resultBlock.classList.add('is-visible');
    });

    let fwd, bkw, prob;

    // move forward
    const Alpha = (state, time, transition_probability, emit_probability, observations, output, fwd) => {
        if (fwd[state][time] !== -1) {
            return fwd[state][time];
        }

        output.push("------------------------------");
        output.push(`Calculating Alpha[State = ${state + 1}, Time = ${time + 1}]`);
        output.push(`Alpha[State, Time] = emit_probability[State][Observation[Time]] * (sum for states from s=1 to ${StatesN} : (Alpha[s, Time - 1]) )`);

        let f = 0;
        let line = "(";
        for(let j = 0; j < StatesN; j++) {
            let add =
            Alpha(j, time - 1, transition_probability, emit_probability, observations, output, fwd) *
                transition_probability[j][state];
            line += add.toPrecision(5);
            if (j !== StatesN - 1) {
                line +=" + ";
            }
            f += add;
        }

        line += ")";

        let mult = emit_probability[state][observations[time]];
        f *= emit_probability[state][observations[time]];

        line = `Alpha[State = ${state + 1}, Time = ${time + 1}] = ${mult.toPrecision(5)} * ${line} = ${f.toPrecision(5)}`;
        output.push(line);
        output.push("------------------------------");

        fwd[state][time] = f;
        return fwd[state][time];
    };

    // move backward
    const Beta = (state, time, transition_probability, emit_probability, observations, output, bkw) => {
        if (bkw[state][time] !== -1) {
            return bkw[state][time];
        }

        output.push("------------------------------");
        output.push(`Calculating Beta[State = ${state + 1}, Time = ${time + 1}]`);
        output.push(`Beta[State, Time] = sum for states from s=1 to ${StatesN} : (Beta[s, Time + 1]) * transition_probability[State][s] * emit_probability[s][Observation[Time + 1]`);

        let b = 0;
        let line = "";
        for(let j = 0; j < StatesN; j++) {

            let add =  Beta(j, time + 1, transition_probability, emit_probability, observations, output, bkw) *
                transition_probability[state][j] * emit_probability[j][observations[time + 1]];
            b += add;
            line += add.toPrecision(5);
            if (j !== StatesN - 1) {
                line +=" + ";
            }
        }

        line = `Beta[State = ${state + 1}, Time = ${time + 1}] = ${line} = ${b.toPrecision(5)}`;

        output.push(line);
        output.push("------------------------------");

        bkw[state][time] = b;
        return bkw[state][time];
    };


    const Forward_Backward = (Y, pi, A, B) => {
        const T = Y.length;
        const K = A.length;

        fwd = getTwoDimensionalArray(StatesN, T, -1);
        bkw = getTwoDimensionalArray(StatesN, T, -1);
        prob = getTwoDimensionalArray(StatesN, T, -1);


        const output = [];
        output.push("------------------------------");
        output.push("Initial computations");
        output.push("For all states from s = 1 to " + StatesN + ":");
        output.push("-- Alpha[s, 1] = emit_probability[s][observations[1]] * П[s]");
        output.push(`-- Beta [s, T] = 1`);
        output.push("                              ");

        for (let i = 0; i < K; i++) {
            fwd[i][0] = B[i][Y[0]] * pi[i];
            bkw[i][T - 1] = 1;
            output.push(`Alpha[${i + 1}][1] = ` + B[i][Y[0]].toPrecision(5) + " * " + pi[i].toPrecision(5) + ' = ' + fwd[i][0].toPrecision(5));
            output.push(`Beta[${i + 1}][${T}]= ` + bkw[i][0].toPrecision(5));
            output.push("------------------------------");

        }
        output.push("Calculate Alpha (Forward)");
        output.push("------------------------------");

        for (let t = 1; t < T; t++) {
            for (let i = 0; i < K; i++) {
                Alpha(i, t, A, B, Y, output, fwd);
            }
        }

        output.push("Calculate Beta (Backward)");
        output.push("------------------------------");

        for (let t = T - 2; t >= 0; t--) {
            for (let i = 0; i < K; i++) {
                Beta(i, t, A, B, Y, output, bkw);
            }
        }

        output.push("Calculate Probabilities* of each State and Time");
        output.push("------------------------------");
        output.push("Prob[State, Time] = Alpha[State][Time] * Beta[State][Time]");
        output.push("We don't normalize our probabilities, so it's ok that they will not equal to 1 in total");

        for (let t = 0; t < T; t++) {
            for (let i = 0; i < K; i++) {
                prob[i][t] = fwd[i][t] * bkw[i][t];
                output.push(`Probability[${i + 1}][${t + 1}]= ` + prob[i][t].toPrecision(5));
                output.push("------------------------------");
            }
        }

        output.push("------------------------------");
        output.push(" Final Sequence");
        let sequence = "[";
        for (let t = 0; t < T; t++) {
            let state = 0;
            for (let i = 1; i < K; i++) {
                if (prob[state][t] < prob[i][t]) {
                    state = i;
                }
            }
            sequence += " " + (state + 1) + " ";
        }

        sequence += "]";
        output.push(sequence);

        return output.join('\n');
    };

    const getObservations = (field) => {
        const chars = field.value.trim().toUpperCase().split('');

        return chars.map((char) => {
            return char.charCodeAt(0) - 65;
        });
    };

    const getInputValues = (element) => {
        const result = [];

        [].forEach.call(element.querySelectorAll('input'), (input) => {
            result.push(parseFloat(input.value));
        });

        return result;
    };

    const getMatrixValues = (element) => {
        const matrix = [];

        [].forEach.call(element.querySelectorAll('.matrix-row'), (row) => {
            matrix.push(getInputValues(row));
        });

        return matrix;
    };

    const generateInputs = (count) => {
        return Array(count)
            .fill('<input class="form-control col-sm-3" type="number" step="any" min="0" max="1" required>')
            .join('\n');
    };

    const generateMatrix = (rows, columns) => {
        const rowHtml = [
            '<div class="row matrix-row">',
                '<li></li>',
                generateInputs(columns),
            '</div>'
        ].join('\n');

        const matrixHtml = [
            '<ol>',
                Array(rows).fill(rowHtml).join('\n'),
            '</ol>'
        ].join('\n');

        return matrixHtml;
    };

    const getTwoDimensionalArray = (rows, columns, value) => {
        return Array(rows).fill(0).map(x => Array(columns).fill(value));
        return Array(rows).fill(Array(columns).fill(value));
    };
});
