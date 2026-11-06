// Filter options
function createAAContainer(containerId) {
    const container = document.getElementById(containerId);
    if (!container) {
        return;
    }
    container.className = 'label-sub';
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];

    container.innerHTML = '';


    for (let i = 1; i <= seqLength; i++) {
        const wrapper = document.createElement('div');
        wrapper.style.flex = '0 0 60px';

        const label = document.createElement('label');
        label.style.color = '#FFF';
        label.style.fontSize = '18px';
        label.style.marginLeft = '1px';
        label.style.cursor = 'pointer';
        label.appendChild(document.createTextNode(`R${i}`));

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = `filterPos`;
        checkbox.value = `R${i}`;
        checkbox.style.cursor = 'pointer';

        label.prepend(checkbox);
        wrapper.appendChild(label);

        const aaGroup = document.createElement('div');
        aaGroup.style.display = 'none';
        aaGroup.style.flexWrap = 'wrap';
        aaGroup.style.marginLeft = '18px';


        aminoAcids.forEach(aa => {
            const aaLabel = document.createElement('label');
            aaLabel.style.color = 'white';
            aaLabel.style.fontSize = '14px';
            aaLabel.style.display = 'flex';
            aaLabel.style.flexDirection = 'column';
            aaLabel.style.alignItems = 'center';
            aaLabel.style.textAlign = 'center';
            aaLabel.style.cursor = 'pointer';

            const aaCheckbox = document.createElement('input');
            aaCheckbox.type = 'checkbox';
            aaCheckbox.name = `${containerId === 'exclAAContainer' ? 'exclR' : 'fixR'}${i}`;
            aaCheckbox.value = aa;
            aaCheckbox.style.cursor = 'pointer';

            aaLabel.appendChild(document.createTextNode(aa));
            aaLabel.appendChild(aaCheckbox);
            aaGroup.appendChild(aaLabel);
        });

        checkbox.addEventListener('change', () => {
            aaGroup.style.display = checkbox.checked ? 'flex' : 'none';
        });

        wrapper.appendChild(aaGroup);
        container.appendChild(wrapper);
    }
}
document.addEventListener('DOMContentLoaded', () => {
    // Monitor AA selection boxes for changes in seqLength
    createAAContainer('fixAAContainer');
    createAAContainer('exclAAContainer');
    createDropContainer();
});
function updateFixedAA() {
    createAAContainer('fixAAContainer');
    createAAContainer('exclAAContainer');
    createDropContainer();
}

function createDropContainer() {
    const container = document.getElementById('dropAAContainer');
    if (!container) return;

    container.className = 'label-sub';
    container.style.display = 'flex';
    container.style.flexWrap = 'wrap';
    container.style.gap = '5px';
    container.style.margin = '5px';

    const seqLength = parseInt(document.getElementById('seqLength').value);
    container.innerHTML = '';

    for (let i = 1; i <= seqLength; i++) {
        const wrapper = document.createElement('div');
        wrapper.style.display = 'flex';
        // wrapper.style.marginTop = '2px';

        const label = document.createElement('label');
        label.style.color = '#FFF';
        label.style.fontSize = '18px';
        label.style.display = 'inline-flex';
        label.style.flexDirection = 'column'; // Stack children vertically
        label.style.cursor = 'pointer';
        label.style.gap = '0px';  // Space between text and checkbox

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = 'dropPos';
        checkbox.value = `R${i}`;

        // Order matters: Text first, then checkbox
        const textNode = document.createTextNode(`R${i}`);
        label.appendChild(textNode);
        label.appendChild(checkbox);

        wrapper.appendChild(label);
        container.appendChild(wrapper);
    }
}



// Combine substrate profiles
function createProfileContainer(subs=true) {
    const container = document.getElementById("profileContainer");
    if (!container) return;

    const nProfiles = parseInt(document.getElementById('nProfiles').value);
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const motifLength = parseInt(document.getElementById('motifLength').value);
    container.innerHTML = '';

    // Header label
    if (subs) {
        container.innerHTML = `
            <div class="help">
                <label>Experimental Data:</label>
                <div class="help-icon">?
                    <span class="help-tooltip">
                        Substrate profiles obtained from "Filter Motif".<br><br>
                        AA Counts:<br>
                        - Count matrix of the substrate profile.<br>
                        - Acceptable file extension: .csv<br><br>
                        Substrates:<br>
                        - Optional input used for a wordcloud.<br>
                        - Acceptable file extension: .json<br><br>
                        Motif Index:<br>
                        - Index of the first AA in the motif within the full substrate sequence.
                    </span>
                </div>
            </div>
        `;
    } else {
        container.innerHTML = `
            <div class="help">
                <label>Experimental Data:</label>
                <div class="help-icon">?
                    <span class="help-tooltip">
                        Substrate profiles obtained from "Filter Motif".<br><br>
                        AA Counts:<br>
                        - Count matrix of the substrate profile.<br>
                        - Acceptable file extension: .csv<br><br>
                        Motif Index:<br>
                        - Index of the first AA in the motif within the full substrate sequence.
                    </span>
                </div>
            </div>
        `;
    }

    const maxIdx = seqLength - motifLength + 1
    for (let i = 1; i <= nProfiles; i++) {
        if (i > maxIdx) {
            break; // Enforce data boundaries
        }
        const l = "121px";

        // Upload file
        const wrapper0 = document.createElement('div');
        wrapper0.innerHTML = `
            <a class="label-w">Set ${i}:</a>
        `;
        container.appendChild(wrapper0);

        const wrapper1 = document.createElement('div');
        wrapper1.className = 'form-wrapper';
        wrapper1.innerHTML = `
            <label class="label-w" style="width: ${l}" for="fileExpCounts${i}">* AA Counts:</label>
            <input type="file" id="fileExpCounts${i}" name="fileExpCounts${i}" accept=".csv" required>
        `;
        container.appendChild(wrapper1);

        if (subs) {
            const wrapper2 = document.createElement('div');
            wrapper2.className = 'form-wrapper';
            wrapper2.innerHTML = `
                <label class="label-w" style="width: ${l};" for="fileExp${i}">* Substrates:</label>
                <input type="file" "id="fileExpSubs${i}" name="fileExp${i}" accept=".json">
            `;
            container.appendChild(wrapper2);
        }

        const wrapper3 = document.createElement('div');
        wrapper3.className = 'form-wrapper';
        wrapper3.style.marginBottom = '1px';
        wrapper3.innerHTML = `
            <label class="label-w" style="width: ${l}" for="fileExp${i}">* Motif Index:</label>
            <input type="number" id="idxStart${i}" name="idxStart${i}" value="${i}"
                   min="1" max="${maxIdx}" style="width: 60px; height: 25px;"
                   required>
        `;
        container.appendChild(wrapper3);
    }
}

document.addEventListener('DOMContentLoaded', () => {
    const page = document.body.getAttribute('data-page');
    if (page === 'combine') {
        createProfileContainer(true);  // Show Substrates
    } else if (page === 'predict') {
        createProfileContainer(false); // Hide Substrates
    }
});

function updateNumProfiles() {
    createProfileContainer(true);
}

function updateNumProfilesPred() {
    checkSubs()
    createProfileContainer(false);
}


async function download() {
    const csrfToken = document.querySelector('meta[name="csrf_token"]').getAttribute('content');
    const response = await fetch('/download', {
        body: JSON.stringify({}),
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    });

    const blob = await response.blob();

    // Extract filename from Content-Disposition header
    let filename = 'comet.zip';
    const disposition = response.headers.get('Content-Disposition');
    if (disposition) {
        const match = disposition.match(/filename[^;=\n]*=([^;\n]*)/);
        if (match && match[1]) {
            filename = match[1].trim().replace(/['"]/g, '');
        }
    }
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}


async function downloadFiles() {
    const csrfToken = document.querySelector('meta[name="csrf_token"]').getAttribute('content');
    const response = await fetch('/downloadFiles', {
        body: JSON.stringify({}),
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    });

    const blob = await response.blob();

    // Extract filename from Content-Disposition header
    let filename = 'TemplateData.zip';
    const disposition = response.headers.get('Content-Disposition');
    console.log('Disposition: ' + disposition)
    if (disposition) {
        const match = disposition.match(/filename[^;=\n]*=([^;\n]*)/);
        if (match && match[1]) {
            filename = match[1].trim().replace(/['"]/g, '');
        }
    }
    console.log('Filename: ' + filename)
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}

async function processForm(formData) {
    const json = {}; // Dont send files as a JSON
    const selectedFixPositions = [];

    // Process the input form
    for (const [key, value] of formData.entries()) {
        if (key === 'filterPos') {
            selectedFixPositions.push(value);  // e.g., ['R2']
        }
        if (json[key]) {
            // When you have more that one value or a key, put the values in a list
            if (!Array.isArray(json[key])) {
                json[key] = [json[key]]; // Convert to array
            }
            json[key].push(value); // Push another value into the list
        } else {
            json[key] = value;
        }
    }

    // Clean out fixR* keys not selected
    Object.keys(json).forEach(key => {
        if (key.startsWith('fix') && !selectedFixPositions.includes(key.replace('fix', ''))) {
            delete json[key];
        }
    });

    // Evaluate job ID
    let jobID = '';
    for (let [key, value] of formData.entries()) {
        if (value instanceof File) {
            if (value.name) {
                jobID += key + '_' + value.name + ' ';
            }
        } else if (value) {
            jobID += key + '_' + value + ' ';
        }
    }
    const date = new Date().toUTCString();
    jobID += date;
    const encoder = new TextEncoder();
    const data = encoder.encode(jobID);
    const hashBuffer = await crypto.subtle.digest('SHA-256', data);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    jobID = hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
    return jobID;
}


//
function pageHome() {
    window.location.href = "/"
}


function pageProcessDNA() {
    window.location.href = "/processDNA";
}


function pageFilterAA() {
    window.location.href = "/filterAA";
}


function pageFilterMotif() {
    window.location.href = "/filterMotif";
}

function pageCombineProfiles() {
    window.location.href = "/combineProfiles";
}

function pagePredictions() {
    window.location.href = "/predictions";
}

// Define button functions
async function buttonProcessDNA() {
    const button = document.querySelector('button[onclick="buttonProcessDNA()"]');
    button.disabled = true;
    const originalText = button.textContent;
    button.textContent = 'Processing';

    const form = document.getElementById("formDNA");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');

    // Evaluate the form
    jobID = await processForm(formData);
    formData.append('jobID', jobID);

    const controller = new AbortController();
    const timeoutMs = 120000;
    console.log('Timeout: ', timeoutMs / 1000, 'sec');
    const timeoutId = setTimeout(() => controller.abort(), timeoutMs);

    try {
        const response = await fetch('/evalFormDNA', {
            body: formData,
            headers: { 'X-CSRFToken': csrfToken },
            method: 'POST',
            signal: controller.signal
        });

        clearTimeout(timeoutId);

        if (response.ok) {
            console.log('Form sent successfully');
            window.location.href = '/results';
        } else {
            const errorText = await response.text();
            console.error(`HTTP ${response.status}:`, errorText);
            alert(`ERROR: ${response.status}.\nCheck console.`);
            button.disabled = false;
            button.textContent = originalText;
        }
    } catch (error) {
        clearTimeout(timeoutId);
        console.error('Fetch Error:', error);
        alert('ERROR: Failed to upload files.\n' +
            'It is recommended to run separate jobs with fewer files.');
        button.disabled = false;
        button.textContent = originalText;
    }
}


async function buttonFilterSubs(filter) {
    // Disable to prevent double click
    const button = document.querySelector('.button');
    button.disabled = true;
    button.textContent = 'Processing';

    // Process form
    const form = document.getElementById("filterSubs");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');

    // Evaluate the form
    if (filter !== 'comet') {
        jobID = await processForm(formData);
        formData.append('jobID', jobID);
    }

    // POST the raw formData to Flask
    if (filter === 'aa') {
        fetch('/evalFormFilterAA', {
            body: formData,  // Send the actual FormData object, not a JSON
            headers: { 'X-CSRFToken': csrfToken },
            method: 'POST'
        })
        .then(response => {
            if (response.ok) {
                window.location.href = '/results'; // Redirect
            } else {
                console.log("ERROR: Filtering substrates.");
                alert("ERROR: Filtering substrates.");
            }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    } else if (filter === 'motif') {
        fetch('/evalFormFilterMotif', {
            body: formData,  // Send the actual FormData object, not a JSON
            headers: { 'X-CSRFToken': csrfToken },
            method: 'POST'
        })
        .then(response => {
            if (response.ok) {
                window.location.href = '/setEntropy'; // Redirect
            } else {
                console.log("ERROR: Filtering motif.");
                alert("ERROR: Filtering motif.");
                }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    } else if (filter === 'comet') {
        fetch('/comet', {
            body: formData,  // Send the actual FormData object, not a JSON
            headers: { 'X-CSRFToken': csrfToken },
            method: 'POST'
        })
        .then(response => {
            if (response.ok) {
                window.location.href = '/results'; // Redirect
            } else {
                console.log("ERROR: Running COMET.");
                alert("ERROR: Running COMET.");
                }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    } else {
        alert('No valid filter type was given to: buttonFilterSubs(filter)')
    }
}

async function buttonCombineProfiles() {
    // Disable to prevent double click
    const button = document.querySelector('.button');
    button.disabled = true;
    button.textContent = 'Processing';

    // Process form
    const form = document.getElementById("filterSubs");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');

    // Evaluate the form
    jobID = await processForm(formData); //
    formData.append('jobID', jobID);

    // POST the raw formData to Flask
    fetch('/evalFormCombineProfiles', {
        body: formData,  // Send the actual FormData object, not a JSON
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    })
    .then(response => {
        if (response.ok) {
            window.location.href = '/results'; // Redirect
        } else {
            console.log("ERROR: Combining profiles.");
            alert("ERROR: Combining profiles.");
        }
    })
    .catch(error => {
        console.error('Error:', error);
        alert("An error occurred.");
    });
}

async function buttonPredActivity() {
    // Disable to prevent double click
    const button = document.querySelector('.button');
    button.disabled = true;
    button.textContent = 'Processing';

    // Process form
    const form = document.getElementById("filterSubs");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');

    // Evaluate the form
    jobID = await processForm(formData); //
    formData.append('jobID', jobID);

    // POST the raw formData to Flask
    fetch('/evalFormPredActivity', {
        body: formData,  // Send the actual FormData object, not a JSON
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    })
    .then(response => {
        if (response.ok) {
            window.location.href = '/results'; // Redirect
        } else {
            console.log("ERROR: Combining profiles.");
            alert("ERROR: Combining profiles.");
        }
    })
    .catch(error => {
        console.error('Error:', error);
        alert("An error occurred.");
    });
}


function addFigure(container, label, fig, fig2 = null) {
    const p = document.createElement('p');
    p.className = 'p2';
    p.textContent = label;
    container.appendChild(p);

    // Add figure
    const img1 = document.createElement('img');
    img1.src = fig;
    img1.style.width = '800px';
    img1.style.maxWidth = '80%';
    img1.style.height = 'auto';
    img1.style.marginBottom = '20px';
    container.appendChild(img1);

    // Add a second figure
    if (fig2) {
        const img2 = document.createElement('img');
        img2.src = fig2;
        img2.style.width = '800px';
        img2.style.maxWidth = '80%';
        img2.style.height = 'auto';
        img2.style.marginBottom = '20px';
        container.appendChild(img2);
    }
}


// Get figures
function getFigures(pollFigs=true, pollInterval) {
    const csrfToken = document.querySelector('meta[name="csrf_token"]').getAttribute('content');
    const container = document.getElementById("figures-container");
    if (!container) return;

    if (pollFigs) {
        const interval = setInterval(() => {
            // new Flask route returning JSON with filenames
            fetch('/checkFigures', {
                method: 'GET',
                headers: { 'X-CSRFToken': csrfToken },
                credentials: 'same-origin'
            })
                .then(res => res.json())
                // Verify if one figure is ready
                .then(data => {
                    container.innerHTML = ''; // Clear figures
                    if (data.entropyProfile) {
                        // Figures: Substrate profile
                        if (data.entropyProfile) addFigure(container, "Substrate Profile: Entropy", data.entropyProfile);
                        if (data.eMapProfile) addFigure(container, "Substrate Profile: Enrichment Map", data.eMapProfile);
                        if (data.eMapScProfile) addFigure(container, "Substrate Profile: Scaled Enrichment Map", data.eMapScProfile);
                        if (data.eLogoProfile && data.eLogoMinProfile) {
                            addFigure(container, "Substrate Profile: Enrichment Logo", data.eLogoProfile, data.eLogoMinProfile);
                        } else {
                            if (data.eLogoProfile) {
                                addFigure(container, "Substrate Profile: Enrichment Logo", data.eLogoProfile);
                            }
                            if (data.eLogoMinProfile) {
                                addFigure(container, "Substrate Profile: Enrichment Logo", data.eLogoMinProfile);
                            }
                        }
                        if (data.wLogoProfile) addFigure(container, "Substrate Profile: WebLogo", data.wLogoProfile);
                        if (data.wordsProfile) addFigure(container, "Substrate Profile: Word Cloud", data.wordsProfile);
                        if (data.barCountsProfile) addFigure(container, "Substrate Profile: Sequence Counts", data.barCountsProfile);
                        if (data.barRFProfile) addFigure(container, "Substrate Profile: Substrate Frequency", data.barRFProfile);
                    } else if (data.predMatrix) {
                        if (data.predMatrix) addFigure(container, 'Prediction Matrix', data.predMatrix);
                        if (data.scatterActivity && data.barPred) {
                            addFigure(container, "Substrate Activity", data.scatterActivity, data.barPred);
                        } else {
                            if (data.scatterActivity) {
                                addFigure(container, "Substrate Activity", data.scatterActivity);
                            }
                            if (data.barPred) {
                                addFigure(container, "Substrate Activity", data.barPred);
                            }
                        }
                        if (data.exp_counts) addFigure(container, "Experimental Counts", data.exp_counts);
                        if (data.bg_counts) addFigure(container, "Background Counts", data.bg_counts);
                    } else if (data.entropy || data.eMap || data.eMapSc ||
                        data.eLogo || data.eLogoMin ||  data.wLogo || data.words ||
                        data.barCounts || data.barRF || data.exp_counts || data.bg_counts) {
                        // Figures
                        if (data.exp_counts) addFigure(container, "Experimental Counts", data.exp_counts);
                        if (data.bg_counts) addFigure(container, "Background Counts", data.bg_counts);
                        if (data.entropy) addFigure(container, 'Entropy', data.entropy);
                        if (data.eMap) addFigure(container, "Enrichment Map", data.eMap);
                        if (data.stDev) addFigure(container, "Standard Deviation", data.stDev);
                        if (data.eMapSc) addFigure(container, "Scaled Enrichment Map", data.eMapSc);
                        if (data.eLogo && data.eLogoMin) {
                            addFigure(container, "Enrichment Logo", data.eLogo, data.eLogoMin);
                        } else {
                            if (data.eLogo) {
                                addFigure(container, "Enrichment Logo", data.eLogo);
                            }
                            if (data.eLogoMin) {
                                addFigure(container, "Enrichment Logo", data.eLogoMin);
                            }
                        }
                        if (data.wLogo) addFigure(container, "WebLogo", data.wLogo);
                        if (data.words) addFigure(container, "Word Cloud", data.words);
                        if (data.barCounts && data.barCountsAll) {
                            addFigure(container, "Sequence Counts", data.barCounts, data.barCountsAll);
                        } else {
                            if (data.barCounts) {
                                addFigure(container, "Sequence Counts", data.barCounts);
                            }
                            if (data.barCountsAll) {
                                addFigure(container, "Sequence Counts", data.barCountsAll);
                            }
                        }
                        if (data.barRF) addFigure(container, "Sequence Frequency", data.barRF);
                    }
                });

            // Only stop polling when job is done
            fetch('/jobStatus')
                .then(r => r.json())
                .then(status => {
                    if (status.jobStatus) {
                        clearInterval(interval);

                        // Append download button to the box
                        const containerBtn = document.getElementById("button-container");
                        if (containerBtn) {
                            const box = document.querySelector('.box');
                            const buttonWrapper = document.createElement('div');
                            buttonWrapper.className = 'button-wrapper';
                            const button = document.createElement('button');
                            button.className = 'button';
                            button.textContent = 'Download';
                            button.onclick = async function() {
                                // button.disabled = true; ##$$
                                // const csrfToken = document.querySelector('meta[name="csrf_token"]').content;
                                button.textContent = 'Downloading';
                                try {
                                    await download();
                                    button.textContent = 'Download';
                                } catch (err) {
                                    console.error('Download error:', err);
                                    button.textContent = 'Download';
                                }
                            };
                            buttonWrapper.appendChild(button);
                            document.getElementById('button-container').appendChild(buttonWrapper);
                        }
                    }
                });
        }, pollInterval); // poll: 1000 = 1 second
    }
}


function getFiguresEntropy(pollFigs=true) {
    const csrfToken = document.querySelector('meta[name="csrf_token"]').getAttribute('content');
    const container = document.getElementById("figures-container");
    if (!container) return;

    if (pollFigs) {
        let finishA = false;
        let finishB = false;
        let finishC = false;
        let finishD = false;

        const interval = setInterval(() => {
            fetch('/checkFigures', {
                method: 'GET',
                headers: { 'X-CSRFToken': csrfToken },
                credentials: 'same-origin'
            })
                .then(res => res.json())
                .then(data => {
                    if (data.entropy && !finishA) {
                        addFigure(container, 'Entropy', data.entropy);
                        finishA = true;
                    }
                    if (data.eMap && !finishB) {
                        addFigure(container, "Enrichment Map", data.eMap);
                        finishB = true;
                    }
                    if (data.eLogo && data.eLogoMin && !finishC && !finishD) {
                        addFigure(container, "Enrichment Logo", data.eLogo, data.eLogoMin);
                        finishC = true;
                        finishD = true;
                    } else {
                        if (data.eLogoProfile && !finishC) {
                            addFigure(container, "Enrichment Logo", data.eLogo);
                            finishC = true;
                        }
                        if (data.eLogoMinProfile && !finishD) {
                            addFigure(container, "Enrichment Logo", data.eLogoMin);
                            finishD = true;
                        }
                    }

                    // Only proceed when BOTH are ready
                    if (finishA && finishB && finishC && finishD) {
                        clearInterval(interval);
                        const minSInput = document.getElementById('minS');
                        const minS = parseFloat(minSInput.value);

                        fetch('/updateMinS', {
                            method: 'POST',
                            headers: {
                                'Content-Type': 'application/json',
                                'X-CSRFToken': csrfToken
                            },
                            body: JSON.stringify({ minS: minS })
                        })
                            .then(res => res.json())
                            .then(updatedData => {
                                updateMotifPosDisplay(updatedData);

                                // Check job status AFTER updateMinS completes
                                return fetch('/jobStatus')
                                    .then(r => r.json());
                            })
                            .then(status => {
                                if (status.jobStatus) {
                                    const buttonContainer = document.getElementById("button-container");
                                    if (buttonContainer) {
                                            const buttonWrapper = document.createElement('div');
                                            buttonWrapper.className = 'button-wrapper';
                                            const button = document.createElement('button');
                                            button.className = 'button button-filter';
                                            button.textContent = 'Filter';
                                            button.onclick = () => buttonFilterSubs('comet');
                                            buttonWrapper.appendChild(button);
                                            buttonContainer.appendChild(buttonWrapper);
                                        }
                                }
                            })
                            .catch(err => console.error('Error in minS update chain:', err));
                    }
                })
                .catch(err => console.error('Error fetching figures:', err));
        }, 1000);
    }
}


function updateMotifPosDisplay(data) {
    const container = document.getElementById('motifPosContainer');
    if (!container) return;
    // remove old entries
    container.querySelectorAll('.motif-pos-entry').forEach(el => el.remove());
    const motifPos = data.motifPos;
    if (motifPos.length > 0) {
        const def = document.getElementById('motifPosDefault');
        if (def) def.style.display = 'none';
        motifPos.forEach(([pos, val]) => {
            const p = document.createElement('p');
            p.className = 'p3 motif-pos-entry';
            if (pos) {
                p.innerHTML = `${pos}: ∆S=<span class="param-value">${val.toFixed(2)}</span>`;
            } else {
                p.innerHTML = `${val}`;
            }
            container.appendChild(p);
        });
    }
}


function updateMinS() {
    const minS = parseFloat(document.getElementById('minS').value);
    const csrfToken = document.querySelector('meta[name="csrf_token"]').getAttribute('content');
    fetch('/updateMinS', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', 'X-CSRFToken': csrfToken },
        body: JSON.stringify({
            minS: minS
        })
    })
    .then(r => r.json())
    .then(data => {
        updateMotifPosDisplay(data);
    });
}


function pollMotifPos() {
    fetch('/motifPos')
        .then(r => r.json())
        .then(data => updateMotifPosDisplay(data));
}


function checkSubs() {
    const seqLength = parseInt(document.getElementById('motifLength').value);
    const rawSubs = document.getElementById(
        'predSubs').value.split(',')
        .map(s => s.trim())
        .filter(s => s !== "");
    const subs = rawSubs.map(item => {
        // If ':' exists, take the part before it; otherwise, take the whole string
        return item.includes(':') ? item.split(':')[0].trim() : item;
    });

    // Inspect inputs
    const errorDiv = document.getElementById('error-message');
    errorDiv.style.display = 'none'; // Hide previous errors
    errorDiv.innerHTML = ''; // Clear content
    const errors = [];

    // Valid amino acids
    const AA = new Set(
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
                'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l',
                'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'
        ]
    );

    // Inspect sequences
    subs.forEach((sub, index) => {
        // 1. Declare notAA here with const so it resets for every sequence
        const notAA = [];

        const isValid = sub.split('').every(char => AA.has(char));

        if (!isValid) {
            sub.split('').forEach((aa, i) => {
                // 2. Add braces for safety and clarity
                if (!AA.has(aa)) {
                    notAA.push(aa);
                }
            });

            // 3. Generate error message
            if (notAA.length === 1) {
                errors.push(`* Sequence #${index + 1}: Contains invalid character "${notAA}"`);
            } else {
                errors.push(`* Sequence #${index + 1}: Contains invalid characters "${notAA.join(", ")}"`);
            }
        }
    });


    subs.forEach((sub, index) => {
         if (sub.length !== seqLength) {
            errors.push(`* Sequence #${index + 1}: ${sub} contains ${sub.length} AAs (expected ${seqLength} AAs)`);
        }
    });

    if (errors.length > 0) {
        errorDiv.innerHTML = "<strong>Validation Failed:</strong><br>" + errors.join("<br>");
        errorDiv.style.display = 'block';
        return false;
    }
    return true;
}
