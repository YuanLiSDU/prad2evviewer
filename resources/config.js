// Reset all frontend state (used by Clear All and mode switching)
function clearFrontend(){
    occData={}; occTcutData={}; occTotal=0;
    eventChannels={}; currentWaveform=null; currentHist={};

    // reset waveform stacking state
    wfStackTraces=[]; wfStackModKey=''; wfStackEnabled=false;
    wfDaqEnabled=false;
    wfRequestId++;  // invalidate any in-flight waveform fetches
    lastHistModule='';
    document.getElementById('wf-stack').checked=false;
    document.getElementById('wf-stack-count').style.display='none';
    document.getElementById('btn-wf-stack-reset').style.display='none';
    document.getElementById('wf-daq').checked=false;
    document.getElementById('wf-daq-info').style.display='none';
    document.getElementById('peaks-table-soft').style.display='';
    document.getElementById('peaks-table-daq').style.display='none';

    // blank DQ plots but keep selected module
    Plotly.react('waveform-div',[], wfLayout(selectedModule?selectedModule.n:'', wfWindowNs()), PC2);
    Plotly.react('heighthist-div',[],{...PL,title:{text:'Height Histogram',font:{size:10,color:'#555'}}},PC2);
    Plotly.react('inthist-div',[],{...PL,title:{text:'Integral Histogram',font:{size:10,color:'#555'}}},PC2);
    Plotly.react('poshist-div',[],{...PL,title:{text:'Peak Position',font:{size:10,color:'#555'}}},PC2);
    document.getElementById('peaks-tbody').innerHTML='';
    document.getElementById('peaks-tbody-daq').innerHTML='';

    // cluster tab
    initClHist(); plotClHist(); plotClStatHists();
    clusterData=null; clusterEvent=-1; selectedCluster=-1;
    currentNclustHist=null; currentNblocksHist=null;
    document.getElementById('cl-select').innerHTML='<option value="all">All</option>';
    document.getElementById('cl-detail-header').innerHTML=
        '<span class="cl-info-text">Click a module or select a cluster</span>';
    document.getElementById('cl-tbody').innerHTML='';

    // LMS tab
    lmsSummaryData=null; lmsSelectedModule=-1; currentLmsData=null;
    _lmsHistRaw=null; _lmsHistModName=null;
    Plotly.react('lms-plot',[],{...PL,title:{text:'LMS History',font:{size:10,color:'#555'}}},PC2);
    document.getElementById('lms-detail-header').innerHTML=
        '<span class="cl-info-text">Click a module to view LMS history</span>';
    document.getElementById('lms-tbody').innerHTML='';
    document.getElementById('lms-ref-select').innerHTML='<option value="-1">None</option>';

    document.getElementById('ring-select').innerHTML='';

    // GEM, EPICS, Physics tabs
    gemEffData=null;
    gemOccupancyData=null;
    // Cluster-tab GEM overlay cache — nullify so the next event refetches
    // and redrawGeo() (below) draws the cluster geo without stale dots.
    gemHits=null; gemHitsEvent=-1;
    // GEM APV waveform tab — drop cached payload, clear the canvas registry,
    // and empty the body so previously rendered traces aren't left visible.
    // Also reset gemApvBuiltKey: buildGemApvSections() uses it as an
    // idempotency guard ("same detector layout → skip rebuild"), so without
    // this reset the next fetch with the same layout would early-return and
    // leave the body empty + gemApvCanvases unpopulated, producing a blank
    // GEM APV auto-report screenshot.
    gemApvData=null;
    gemApvCanvases.clear();
    gemApvBuiltKey='';
    gemApvSkippedSincePause=0;
    const apvBody=document.getElementById('gem-apv-body');
    if(apvBody) apvBody.innerHTML='';
    clearEpicsFrontend();
    clearPhysicsFrontend();

    sampleCount=0;
    updateHeaderStats();
    redrawGeo();
    document.getElementById('status-bar').textContent='All data cleared';
}

// fetch /api/config and reconfigure the UI.  Returns the underlying
// promise so callers (e.g. the Cut-Settings Save handler) can await
// applyConfig() before triggering plot redraws that depend on the
// freshly-updated histConfig fields.
function fetchConfigAndApply(){
    return fetch('/api/config').then(r=>r.json()).then(applyConfig);
}

// ET connection dialog
function openEtDialog(){
    // populate fields with current ET config from last /api/config
    const etc=window._etConfig||{};
    document.getElementById('et-input-host').value=etc.host||'localhost';
    document.getElementById('et-input-port').value=etc.port||11111;
    document.getElementById('et-input-file').value=etc.et_file||'/tmp/et_sys_prad2';
    document.getElementById('et-input-station').value=etc.station||'prad2_monitor';
    document.getElementById('et-status-msg').textContent='';
    document.getElementById('et-backdrop').classList.add('open');
    document.getElementById('et-dialog').classList.add('open');
}
function closeEtDialog(){
    document.getElementById('et-backdrop').classList.remove('open');
    document.getElementById('et-dialog').classList.remove('open');
}

function applyConfig(data){
    const crateRoc=data.crate_roc||{};
    const rawMap=data.hycal_map||[];
    modules=[];
    for(const m of rawMap){const g=m.geo,d=m.daq;if(!g||!d)continue;
        modules.push({n:m.n,t:m.t==='PbGlass'?'G':'W',
            x:g.x,y:g.y,sx:g.sx,sy:g.sy,
            roc:crateRoc[String(d.crate)]||0,sl:d.slot,ch:d.channel});}
    totalEvents=data.total_events||0;
    filterActive=data.filter_active||false;
    filteredCount=data.filtered_count||totalEvents;
    histEnabled=data.hist_enabled||false;
    histConfig=data.hist||{};
    // Waveform-Tab peak filter — runtime-mutable via the Cut-Settings dialog.
    // `waveform_filter_active` is the server-side enable flag (the "apply"
    // toggle); `waveform_filter` is the {time, integral, height, quality_bits}
    // payload.  Distinct from `filter_active` (event-level filter loaded via
    // /api/filter/load).
    histConfig.waveform_filter         = data.waveform_filter         || {};
    histConfig.waveform_filter_active  = !!data.waveform_filter_active;
    histConfig.waveform_filter_default = data.waveform_filter_default || {};
    histConfig.quality_bits            = data.quality_bits            || [];
    // sync "apply" toggle to server's enable state
    const cutApplyCb = document.getElementById('cut-apply');
    if (cutApplyCb) cutApplyCb.checked = histConfig.waveform_filter_active;
    refLines=data.ref_lines||{};
    triggerBitsDef=data.trigger_bits||[];
    triggerTypeDef=data.trigger_type||[];
    // load per-tab trigger filters from server config
    const tf=data.trigger_filter||{};
    for(const [tab, filt] of Object.entries(tf)){
        tabTrigFilter[tab]={accept:filt.trigger_accept||0, reject:filt.trigger_reject||0};
    }
    buildTriggerFilterUI();
    restoreTrigFilterFromTab();
    // cluster histogram configs
    if(data.cluster_hist){
        clHistMin=data.cluster_hist.min||0;
        clHistMax=data.cluster_hist.max||3000;
        clHistStep=data.cluster_hist.step||10;
    }
    if(data.nclusters_hist){
        nclustMin=data.nclusters_hist.min||0;
        nclustMax=data.nclusters_hist.max||20;
        nclustStep=data.nclusters_hist.step||1;
    }
    if(data.nblocks_hist){
        nblocksMin=data.nblocks_hist.min||0;
        nblocksMax=data.nblocks_hist.max||40;
        nblocksStep=data.nblocks_hist.step||1;
    }
    if(data.raw_energy_hist){
        rawEnergyMin=data.raw_energy_hist.min||0;
        rawEnergyMax=data.raw_energy_hist.max||6000;
        rawEnergyStep=data.raw_energy_hist.step||20;
    }
    initClHist();
    if(data.lms){
        g_lmsWarnThresh=data.lms.warn_threshold||0.1;
        const sel=document.getElementById('lms-ref-select');
        sel.innerHTML='<option value="-1">None</option>';
        if(data.lms.ref_channels){
            for(const rc of data.lms.ref_channels){
                const o=document.createElement('option');
                o.value=rc.index;
                o.textContent=rc.name;
                sel.appendChild(o);
            }
        }
        // Drift baseline indicator pill — short label shows the run number
        // (extracted from the .dat filename, e.g. "prad_024352_LMS.dat" →
        // "Baseline: run 24352").  Tooltip carries the full path + per-type
        // bands + suppress list so operators can audit what's being applied
        // without leaving the LMS tab.
        const pill=document.getElementById('lms-baseline-pill');
        if(pill){
            if(data.lms.drift_enabled){
                const path=data.lms.drift_baseline||'';
                const base=path.split('/').pop()||path;
                const m=base.match(/_0*(\d+)_/);
                const label=m?`run ${m[1]}`:base;
                pill.style.display='inline-block';
                pill.textContent=`Baseline: ${label}`;
                const fmt=(lo,hi)=>(lo!=null&&hi!=null)?`[${(+lo).toFixed(2)}, ${(+hi).toFixed(2)}]`:'?';
                const sup=(data.lms.drift_suppress_types||[]).join(', ')||'(none)';
                pill.title =
                    `Drift detection baseline\n`+
                    `File:        ${path}\n`+
                    `Ref channel: ${data.lms.drift_ref||'(default)'}\n`+
                    `W band:      ${fmt(data.lms.drift_low_w, data.lms.drift_high_w)}\n`+
                    `G band:      ${fmt(data.lms.drift_low_g, data.lms.drift_high_g)}\n`+
                    `Suppressed:  ${sup}`;
            } else {
                pill.style.display='none';
                pill.textContent='Baseline: —';
                pill.title='Drift detection disabled (no baseline configured)';
            }
        }
    }
    if(data.color_ranges){
        for(const [k,v] of Object.entries(data.color_ranges)){
            if(Array.isArray(v) && v.length===2 && !geoRangeOverrides[k])
                geoRangeOverrides[k]=v;
        }
    }
    if(data.refresh_ms){
        refreshEventMs=data.refresh_ms.event||200;
        refreshRingMs=data.refresh_ms.ring||500;
        refreshHistMs=data.refresh_ms.histogram||2000;
        refreshLmsMs=data.refresh_ms.lms||2000;
    }
    if(data.monitor_status){
        const lt=data.monitor_status.livetime||{};
        const beam=data.monitor_status.beam||{};
        const be=beam.energy||{};
        const bc=beam.current||{};
        livetimeEnabled    = !!(lt.enabled || lt.measured_enabled);
        livetimeUnit       = lt.unit || '%';
        livetimeHealthy    = lt.healthy ?? 90;
        livetimeWarning    = lt.warning ?? 80;
        beamEnergyEnabled  = !!be.enabled;
        beamEnergyUnit     = be.unit || 'MeV';
        beamCurrentEnabled = !!bc.enabled;
        beamCurrentUnit    = bc.unit || 'nA';
        beamCurrentTripWarn= (bc.trip_warn_below!=null) ? bc.trip_warn_below : null;
        // Use the smallest configured poll_sec so the snappier metric drives
        // the loop; floor at 1s to keep the server from being hammered.
        const polls=[lt.poll_sec, be.poll_sec, bc.poll_sec]
            .filter(v=>v!=null && v>0);
        const minPoll=polls.length ? Math.min(...polls) : 5;
        monitorStatusPollMs   = Math.max(1000, minPoll*1000);
        monitorStatusEnabled  = livetimeEnabled || beamEnergyEnabled || beamCurrentEnabled;
    }
    initReport(data);
    initEpics(data);
    initPhysics(data);
    updateTimeCutLabel();
    mode=data.mode||'file';
    etAvailable=data.et_available||false;
    if(data.et_config) window._etConfig=data.et_config;
    fileAvailable=data.file_available||false;
    if(data.source) sourceCaps=data.source;
    const appTitle=mode==='online'?'PRad-II HyCal Monitor':'PRad-II HyCal Event Viewer';
    document.title=appTitle;
    document.getElementById('app-title').textContent=appTitle;
    g_currentFile=data.current_file||'';
    g_dataDirEnabled=data.data_dir_enabled||false;
    g_dataDir=data.data_dir||'';
    g_histCheckbox=histEnabled;

    const hcb=document.getElementById('hist-checkbox');
    if(hcb) hcb.checked=histEnabled;

    // show/hide mode-specific UI
    document.getElementById('nav-file').style.display   = mode!=='online'?'flex':'none';
    document.getElementById('nav-online').style.display = mode==='online'?'flex':'none';
    // Monitor status (livetime + beam) — start/stop polling per mode.
    if(mode==='online') startMonitorStatusPolling(); else stopMonitorStatusPolling();
    document.getElementById('btn-open').style.display='';

    // mode toggle button — visible whenever ET is available
    const toggleBtn=document.getElementById('btn-mode-toggle');
    if(etAvailable){
        toggleBtn.style.display='';
        toggleBtn.textContent=mode==='online'?'Go Offline':'Go Online';
    } else {
        toggleBtn.style.display='none';
    }

    // hide tabs based on data source capabilities
    document.querySelectorAll('.tab').forEach(t=>{
        const tab=t.dataset.tab;
        if(tab==='lms')    t.style.display=sourceCaps.has_waveforms?'':'none';
        if(tab==='epics')  t.style.display=sourceCaps.has_epics?'':'none';
    });

    // auto-switch to cluster tab if source has no waveform data
    if(!sourceCaps.has_waveforms && activeTab==='dq'){
        switchTab('cluster');
    }

    if(mode==='file'){
        if(filterActive){
            document.getElementById('ev-total').textContent=`/ ${filteredCount} (${totalEvents} total)`;
            fetch('/api/filter/indices').then(r=>r.json()).then(idx=>{
                filteredIndices=idx;
                if(idx.length>0) loadEvent(idx[0]);
            }).catch(()=>{filteredIndices=null;});
        } else {
            document.getElementById('ev-total').textContent=`/ ${totalEvents}`;
            filteredIndices=null;
        }

        updateHeaderStats();
        if(histEnabled) { fetchOccupancy(); fetchClHist(); fetchGemResiduals(); }
        fetchEpicsChannels(); fetchEpicsLatest();
        if(activeTab==='epics') fetchAllEpicsSlots();
        if(activeTab==='physics') fetchPhysics();
        syncDqRange();
        geoViewInit=false; resizeGeo();
        if(totalEvents>0 && !filterActive) loadEvent(1);
    } else if(mode==='online'){
        setEtStatus(data.et_connected||false, !data.et_connected);
        syncDqRange();
        fetchOccupancy();
        fetchEpicsChannels(); fetchEpicsLatest();
        if(activeTab==='physics') fetchPhysics();
        resizeGeo();
        updateRingSelector();
        loadLatestEvent();
    } else {
        // idle mode
        syncDqRange();
        resizeGeo();
        updateHeaderStats();
    }
}
