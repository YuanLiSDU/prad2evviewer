// lms.js — LMS gain monitoring tab: geo provider, history, summary table

let g_lmsWarnThresh=0.1;
let g_lmsRefIndex=-1;  // -1 = None (no normalization)
let currentLmsData=null;  // {x:[], y:[]} for copy button
let lmsSummaryData=null;  // {modules:{idx:{name,mean,rms,count,warn}}, events}
let lmsSelectedModule=-1;
// Last /api/lms/<idx> response + module label, kept so theme flips can
// re-render the lms-plot from cache (Plotly bakes THEME.success/danger/text
// /overlay into trace lines and the legend bg at draw time).
let _lmsHistRaw=null, _lmsHistModName=null;

// Tri-state from the API: 'drift' (top-priority error from gain-drift
// detection) > 'warn' (rms/floor stability check) > 'ok'.  Older API
// responses may only have md.warn — fall back to it so the GUI keeps
// working against an old server.
function lmsState(md){
    if(!md) return 'ok';
    if(md.state) return md.state;
    return md.warn ? 'warn' : 'ok';
}

const LMS_DRIFT_SUFFIX_RE = /\s*\(\d+ drift\)$/;
function setPillText(pill, text){
    if(pill.textContent !== text) pill.textContent = text;
}

function geoLms(){
    const metric=document.getElementById('lms-color-metric').value;
    const useLog=document.getElementById('lms-log-scale').checked;
    const mods=lmsSummaryData?lmsSummaryData.modules:{};

    // For 'drift' metric: paint distance from 1.0 so both gain loss (drift<1)
    // and gain growth (drift>1) light up.  Modules with no baseline contribute
    // 0 (= empty colour for that module type).
    const lmsVal=md=>{
        if(!md) return null;
        if(metric==='mean')     return md.mean;
        if(metric==='rms')      return md.rms;
        if(metric==='rms_frac') return md.mean>0?md.rms/md.mean:0;
        if(metric==='drift')    return (md.drift!=null) ? Math.abs(md.drift-1) : 0;
        // 'warn' is the catch-all status colorizer.
        return lmsState(md)==='ok' ? 0 : 1;
    };

    let autoMax=0;
    for(const k in mods){ const v=lmsVal(mods[k]); if(v>autoMax) autoMax=v; }
    if(autoMax<=0) autoMax=1;
    const lmsr=getGeoRange('lms',metric);
    const vmin=lmsr[0]!==null?lmsr[0]:0;
    const vmax=lmsr[1]!==null?lmsr[1]:autoMax;
    const dp=(metric==='rms_frac'||metric==='drift')?3:metric==='rms'?2:0;
    document.getElementById('lms-range-min-show').textContent=vmin.toFixed(dp);
    document.getElementById('lms-range-max-show').textContent=vmax.toFixed(dp);

    renderGeo(
        i => {
            const md=mods[String(i)];
            // Status metric: drift = red, warn = orange, ok = green
            if(metric==='warn'){
                const st=lmsState(md);
                if(st==='drift') return THEME.danger;
                if(st==='warn')  return THEME.warn || THEME.danger;
                if(md && md.count>0) return THEME.success;
                return geoEmptyColor(modules[i].t);
            }
            const val=lmsVal(md);
            if(val!==null&&val>0){
                return geoValueColor(val,vmin,vmax,useLog);
            }
            return geoEmptyColor(modules[i].t);
        },
        i => {
            if(lmsSelectedModule===i) return {color:THEME.selectBorder,width:2.5};
            const md=mods[String(i)];
            const st=lmsState(md);
            if(st==='drift') return {color:THEME.danger,width:2};
            if(st==='warn')  return {color:THEME.warn||THEME.danger,width:1.5};
            return null;
        },
        null
    );
}

function fetchLmsSummary(){
    const refQ=g_lmsRefIndex>=0?`?ref=${g_lmsRefIndex}`:'';
    return fetch(`/api/lms/summary${refQ}`).then(r=>r.json()).then(data=>{
        lmsSummaryData=data;
        geoLms();
        updateLmsTable();
        if(hoveredModule) updateGeoTooltip();
        // Always refresh dot+pill state.  The dot itself is suppressed when
        // on the LMS tab (see below) — but the pill lives in the LMS toolbar
        // and would otherwise freeze with stale drift-count text.
        updateLmsDot();
    }).catch(()=>{});
}

function updateLmsDot(){
    const dot=document.getElementById('lms-dot');
    const pill=document.getElementById('lms-baseline-pill');
    if(!lmsSummaryData||!lmsSummaryData.modules){
        if(activeTab!=='lms') dot.className='tab-dot';
        if(pill){
            pill.classList.remove('has-drift');
            setPillText(pill, pill.textContent.replace(LMS_DRIFT_SUFFIX_RE,''));
        }
        return;
    }
    // Drift = top-priority error → solid red ('alert').  Warn-only =
    // orange dot ('warn').  All-ok = no dot.
    let hasDrift=false, hasWarn=false, driftCount=0;
    for(const m of Object.values(lmsSummaryData.modules)){
        const st=lmsState(m);
        if(st==='drift'){ hasDrift=true; driftCount++; }
        else if(st==='warn') hasWarn=true;
    }
    // Dot is for the *inactive* tab — switchTab clears it on entry, and
    // re-asserting alert/warn while the user is looking at the data would
    // be visual noise.
    if(activeTab!=='lms')
        dot.className='tab-dot'+(hasDrift?' alert':hasWarn?' warn':'');
    // Pill flashes red when drift is outstanding and gets the count appended,
    // so operators see "problem + which baseline" in one glance.
    if(pill && pill.style.display!=='none'){
        pill.classList.toggle('has-drift', hasDrift);
        const base=pill.textContent.replace(LMS_DRIFT_SUFFIX_RE,'');
        setPillText(pill, hasDrift ? `${base} (${driftCount} drift)` : base);
    }
}

function fetchLmsHistory(modIdx, modName){
    const refQ=g_lmsRefIndex>=0?`?ref=${g_lmsRefIndex}`:'';
    fetch(`/api/lms/${modIdx}${refQ}`).then(r=>r.json()).then(data=>{
        _lmsHistRaw=data;
        _lmsHistModName=modName;
        renderLmsHistory();
    }).catch(()=>{});
}

function renderLmsHistory(){
    if(!_lmsHistRaw || _lmsHistModName==null) return;
    const data=_lmsHistRaw;
    const modName=_lmsHistModName;
    if(!data.time||!data.time.length){
        currentLmsData=null;
        Plotly.react('lms-plot',[],{...PL,
            title:{text:`${modName} — No LMS data`,font:{size:10,color:THEME.textMuted}}},PC2);
        return;
    }
    currentLmsData={x:Array.from(data.time), y:Array.from(data.integral)};
    const vals=data.integral;
    const mean=vals.reduce((a,b)=>a+b,0)/vals.length;
    const warnHi=mean*(1+g_lmsWarnThresh);
    const warnLo=mean*(1-g_lmsWarnThresh);
    const tRange=[data.time[0],data.time[data.time.length-1]];

    Plotly.react('lms-plot',[
        {x:data.time, y:data.integral, type:'scatter', mode:'markers',
         marker:{color:'#ff922b',size:3}, name:'LMS integral'},
        {x:tRange, y:[mean,mean],
         type:'scatter', mode:'lines', line:{color:THEME.success,width:1,dash:'dash'}, name:`Mean ${mean.toFixed(0)}`},
        {x:tRange, y:[warnHi,warnHi],
         type:'scatter', mode:'lines', line:{color:THEME.danger,width:1,dash:'dot'}, showlegend:false},
        {x:tRange, y:[warnLo,warnLo],
         type:'scatter', mode:'lines', line:{color:THEME.danger,width:1,dash:'dot'}, showlegend:false},
    ],{...PL,
        title:{text:`LMS — ${modName} (${data.events} pts)${g_lmsRefIndex>=0?' [ref corrected]':''}`,
            font:{size:10,color:THEME.text}},
        xaxis:{...PL.xaxis,
            title:data.sync_unix
                ?`Time (s) after ${new Date((data.sync_unix - data.sync_rel_sec)*1000).toISOString().replace('T',' ').slice(0,19)} UTC`
                :'Time (s)'},
        yaxis:{...PL.yaxis,title:g_lmsRefIndex>=0?'Corrected Integral':'Integral'},
        legend:{x:1,y:1,xanchor:'right',bgcolor:THEME.overlay,font:{size:9}},
        margin:{...PL.margin,t:28,b:36},
        shapes:refShapes('lms'),
    },PC2);

    document.getElementById('lms-info-text').innerHTML=
        `<span class="mod-name">${modName}</span> <span class="mod-daq">Mean: ${mean.toFixed(1)} | RMS: ${(Math.sqrt(vals.reduce((s,v)=>s+(v-mean)**2,0)/vals.length)).toFixed(1)} | ${data.events} pts</span>`;
}

function updateLmsTable(){
    const tbody=document.getElementById('lms-tbody');
    if(!lmsSummaryData||!lmsSummaryData.modules){
        tbody.innerHTML='<tr><td colspan="7" style="text-align:center;color:var(--dim);padding:8px">No LMS data</td></tr>';
        return;
    }
    // Sort: drift first (largest |drift-1| at top), then warn (largest
    // RMS/Mean), then OK.  This puts the most-broken channels on top so
    // the operator sees them without scrolling.
    const entries=Object.entries(lmsSummaryData.modules).map(([idx,m])=>({idx:parseInt(idx),...m}));
    const stateRank=st=>(st==='drift'?0:st==='warn'?1:2);
    entries.sort((a,b)=>{
        const ra=stateRank(lmsState(a)), rb=stateRank(lmsState(b));
        if(ra!==rb) return ra-rb;
        if(lmsState(a)==='drift'){
            // Worst drift first within the drift group.
            const da=a.drift!=null?Math.abs(a.drift-1):0;
            const db=b.drift!=null?Math.abs(b.drift-1):0;
            return db-da;
        }
        const ramf=a.mean>0?a.rms/a.mean:0, rbmf=b.mean>0?b.rms/b.mean:0;
        return rbmf-ramf;
    });
    let rows='';
    for(const e of entries){
        const rmsFrac=e.mean>0?(e.rms/e.mean*100).toFixed(1):'--';
        const sel=lmsSelectedModule===e.idx;
        const st=lmsState(e);
        // Drift cell: show value, mark with a (·) hint when suppressed so the
        // operator knows the module is out of band but the warn was muted on
        // purpose.  No row tint — suppressed rows stay 'ok' in sort order.
        let driftCell;
        if(e.drift==null){
            driftCell='--';
        } else if(e.drift_suppressed){
            driftCell=`<span title="warning suppressed by type" style="color:var(--dim)">${e.drift.toFixed(2)}·</span>`;
        } else {
            driftCell=e.drift.toFixed(2);
        }
        let statusHtml;
        if(st==='drift')      statusHtml='<span class="lms-drift">DRIFT</span>';
        else if(st==='warn')  statusHtml='<span class="lms-warn">WARN</span>';
        else                  statusHtml='<span class="lms-ok">OK</span>';
        const rowCls='cl-table-row'+(sel?' selected':'')+(st==='drift'?' drift-row':'');
        rows+=`<tr class="${rowCls}" data-idx="${e.idx}">
            <td>${e.name}</td>
            <td>${e.mean.toFixed(1)}</td>
            <td>${e.rms.toFixed(2)}</td>
            <td>${rmsFrac}%</td>
            <td style="text-align:center">${driftCell}</td>
            <td style="text-align:center">${e.count}</td>
            <td style="text-align:center">${statusHtml}</td>
        </tr>`;
    }
    tbody.innerHTML=rows;
    tbody.querySelectorAll('.cl-table-row').forEach(tr=>{
        tr.onclick=()=>{
            const idx=parseInt(tr.dataset.idx);
            lmsSelectedModule=idx;
            const mod=modules.find(m=>modules.indexOf(m)===idx);
            const name=lmsSummaryData.modules[idx]?lmsSummaryData.modules[idx].name:'';
            fetchLmsHistory(idx, name);
            updateLmsTable();
            geoLms();
        };
    });
}

// Theme flip — lms-plot embeds THEME.success/danger/text/overlay/textMuted
// in trace lines, the legend bg, and titles.  Replay from the cached raw
// response (set by fetchLmsHistory) so a flip doesn't leave stale colors.
// geoLms paints the SHARED geo canvas — viewer.js's master listener calls
// redrawGeo() for the active tab, so we don't repeat it here.
if (typeof onThemeChange === 'function') {
    onThemeChange(renderLmsHistory);
}
