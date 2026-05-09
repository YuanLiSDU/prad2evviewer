// physics.js — Physics tab: HyCal cluster XY (left) + Møller XY (top right) + energy vs angle (bottom right)
//
// Depends on globals from viewer.js: PL, PC_EPICS, activeTab

let physicsData=null, mollerData=null, hycalXyData=null;

function fetchEnergyAngle(){
    return fetch('/api/physics/energy_angle').then(r=>r.json()).then(data=>{
        physicsData=data;
        plotEnergyAngle();
    }).catch(()=>{});
}

function fetchMoller(){
    return fetch('/api/physics/moller').then(r=>r.json()).then(data=>{
        mollerData=data;
        plotMollerXY();
    }).catch(()=>{});
}

function fetchHycalXY(){
    return fetch('/api/physics/hycal_xy').then(r=>r.json()).then(data=>{
        hycalXyData=data;
        plotHycalXY();
    }).catch(()=>{});
}

function fetchPhysics(){
    return Promise.all([fetchEnergyAngle(), fetchMoller(), fetchHycalXY()]);
}

// ep elastic scattering: E' = E / (1 + (E/Mp)*(1 - cos(theta)))
function elasticEp(beamE, thetaDeg){
    const Mp=938.272;
    const th=thetaDeg*Math.PI/180;
    return beamE/(1+(beamE/Mp)*(1-Math.cos(th)));
}

function plotEnergyAngle(){
    const div='physics-plot';
    if(!physicsData||!physicsData.bins||!physicsData.bins.length||!physicsData.nx){
        Plotly.react(div,[],{...PL,title:{text:'Energy vs Angle — No data',font:{size:12,color:THEME.textDim}}},PC_EPICS);
        document.getElementById('physics-stats').textContent='';
        return;
    }
    const d=physicsData;
    const logZ=document.getElementById('physics-logz').checked;
    const showElastic=document.getElementById('physics-elastic').checked;

    const z=[];
    for(let iy=0;iy<d.ny;iy++){
        const row=d.bins.slice(iy*d.nx,(iy+1)*d.nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
    }
    const x=[];for(let i=0;i<d.nx;i++) x.push(d.angle_min+(i+0.5)*d.angle_step);
    const y=[];for(let i=0;i<d.ny;i++) y.push(d.energy_min+(i+0.5)*d.energy_step);

    const traces=[{
        z:z, x:x, y:y,
        type:'heatmap', colorscale:'Hot', reversescale:false,
        hovertemplate:'θ=%{x:.2f}° E=%{y:.0f} MeV: %{text}<extra></extra>',
        text:z.map((row,iy)=>row.map((v,ix)=>String(d.bins[iy*d.nx+ix]))),
        colorbar:{title:logZ?'log₁₀(counts)':'counts',titleside:'right',
            titlefont:{size:10,color:THEME.textDim},tickfont:{size:9,color:THEME.textDim}},
    }];

    if(showElastic && d.beam_energy>0){
        const ex=[],ey=[];
        for(let th=d.angle_min+0.1;th<=d.angle_max;th+=0.05){
            const e=elasticEp(d.beam_energy,th);
            if(e>=d.energy_min&&e<=d.energy_max){ex.push(th);ey.push(e);}
        }
        traces.push({x:ex,y:ey,mode:'lines',
            line:{color:THEME.success,width:2,dash:'dot'},
            name:`ep elastic (${d.beam_energy.toFixed(2)} MeV)`,
            hovertemplate:'θ=%{x:.2f}° E=%{y:.0f} MeV<extra>ep elastic</extra>'});
    }

    Plotly.react(div,traces,{...PL,
        title:{text:`Energy vs Angle (${d.events} evts)`,font:{size:12,color:THEME.text}},
        xaxis:{...PL.xaxis,title:'Scattering Angle (deg)'},
        yaxis:{...PL.yaxis,title:'Energy (MeV)'},
        margin:{l:55,r:80,t:30,b:40},
        showlegend:showElastic,
        legend:{x:0.7,y:0.95,font:{size:10,color:THEME.textDim},bgcolor:'rgba(0,0,0,0)'},
        shapes:refShapes('energy_angle'),
    },PC_EPICS);

    // stats line
    const ml=mollerData;
    let stats=`${d.events} evts | beam: ${d.beam_energy>0?d.beam_energy.toFixed(2):'?'} MeV`;
    if(ml) stats+=` | Møller: ${ml.moller_events}`;
    document.getElementById('physics-stats').textContent=stats;
}

function plotMollerXY(){
    const div='moller-xy-plot';
    const d=mollerData;
    if(!d||!d.xy_bins||!d.xy_bins.length||!d.xy_nx){
        Plotly.react(div,[],{...PL,title:{text:'Møller XY — No data',font:{size:12,color:THEME.textDim}}},PC_EPICS);
        return;
    }
    const logZ=document.getElementById('physics-logz').checked;
    const z=[];
    for(let iy=0;iy<d.xy_ny;iy++){
        const row=d.xy_bins.slice(iy*d.xy_nx,(iy+1)*d.xy_nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
    }
    const x=[];for(let i=0;i<d.xy_nx;i++) x.push(d.xy_x_min+(i+0.5)*d.xy_x_step);
    const y=[];for(let i=0;i<d.xy_ny;i++) y.push(d.xy_y_min+(i+0.5)*d.xy_y_step);

    const cuts=d.cuts||{};
    const fmtA=v=>v!=null?v.toFixed(2):'?';
    const cutTxt=`θ∈[${fmtA(cuts.angle_min)},${fmtA(cuts.angle_max)}]° Esum±${((cuts.energy_tolerance||0.1)*100).toFixed(0)}%`;

    // θ ring overlay: convert the moller angle window into HyCal-plane radii
    // via r = dz · tan(θ), centered at (target_x, target_y).  dz is the
    // target→HyCal lever arm — same geometry the server uses to compute the
    // per-cluster theta in app_state.cpp.
    const shapes=[...(refShapes('moller_xy')||[])];
    const target=d.target||[0,0,0];
    const dz=((d.hycal_z!=null?d.hycal_z:0)-target[2]);
    if(dz>0 && cuts.angle_min!=null && cuts.angle_max!=null){
        const cx=target[0], cy=target[1];
        const ringColor=THEME.accent;
        [cuts.angle_min,cuts.angle_max].forEach(thDeg=>{
            const r=dz*Math.tan(thDeg*Math.PI/180);
            shapes.push({
                type:'circle', xref:'x', yref:'y',
                x0:cx-r, x1:cx+r, y0:cy-r, y1:cy+r,
                line:{color:ringColor,width:1.2,dash:'dash'},
                fillcolor:'rgba(0,0,0,0)',
            });
        });
    }

    Plotly.react(div,[{
        z:z, x:x, y:y,
        type:'heatmap', colorscale:'Hot', reversescale:false,
        hovertemplate:'x=%{x:.1f} y=%{y:.1f} mm: %{text}<extra></extra>',
        text:z.map((row,iy)=>row.map((v,ix)=>String(d.xy_bins[iy*d.xy_nx+ix]))),
        colorbar:{title:logZ?'log₁₀':'counts',titleside:'right',
            titlefont:{size:10,color:THEME.textDim},tickfont:{size:9,color:THEME.textDim}},
    }],{...PL,
        title:{text:`Møller XY (${d.moller_events} evts) ${cutTxt}`,font:{size:11,color:THEME.text}},
        xaxis:{...PL.xaxis,title:'X (mm)',scaleanchor:'y',scaleratio:1},
        yaxis:{...PL.yaxis,title:'Y (mm)'},
        margin:{l:50,r:70,t:30,b:35},
        shapes:shapes,
    },PC_EPICS);
}

function plotHycalXY(){
    const div='hycal-xy-plot';
    const d=hycalXyData;
    if(!d||!d.xy_bins||!d.xy_bins.length||!d.xy_nx){
        Plotly.react(div,[],{...PL,title:{text:'HyCal Cluster Hits — No data',font:{size:12,color:THEME.textDim}}},PC_EPICS);
        return;
    }
    const logZ=document.getElementById('physics-logz').checked;
    const z=[];
    for(let iy=0;iy<d.xy_ny;iy++){
        const row=d.xy_bins.slice(iy*d.xy_nx,(iy+1)*d.xy_nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
    }
    const x=[];for(let i=0;i<d.xy_nx;i++) x.push(d.xy_x_min+(i+0.5)*d.xy_x_step);
    const y=[];for(let i=0;i<d.xy_ny;i++) y.push(d.xy_y_min+(i+0.5)*d.xy_y_step);

    const c=d.cuts||{};
    const fracPct=((c.energy_frac_min||0.9)*100).toFixed(0);
    const cutTxt=`Ncl=${c.n_clusters||1}, E≥${fracPct}% Eb, blocks∈[${c.nblocks_min||0},${c.nblocks_max||0}]`;

    Plotly.react(div,[{
        z:z, x:x, y:y,
        type:'heatmap', colorscale:'Hot', reversescale:false,
        hovertemplate:'x=%{x:.1f} y=%{y:.1f} mm: %{text}<extra></extra>',
        text:z.map((row,iy)=>row.map((v,ix)=>String(d.xy_bins[iy*d.xy_nx+ix]))),
        colorbar:{title:logZ?'log₁₀':'counts',titleside:'right',
            titlefont:{size:10,color:THEME.textDim},tickfont:{size:9,color:THEME.textDim}},
    }],{...PL,
        title:{text:`HyCal Cluster Hits (${d.events} evts) ${cutTxt}`,font:{size:11,color:THEME.text}},
        xaxis:{...PL.xaxis,title:'X (mm)',scaleanchor:'y',scaleratio:1},
        yaxis:{...PL.yaxis,title:'Y (mm)'},
        margin:{l:50,r:70,t:30,b:35},
        shapes:refShapes('hycal_xy'),
    },PC_EPICS);
}

function clearPhysicsFrontend(){
    physicsData=null; mollerData=null; hycalXyData=null;
    Plotly.react('physics-plot',[],{...PL},PC_EPICS);
    Plotly.react('moller-xy-plot',[],{...PL},PC_EPICS);
    Plotly.react('hycal-xy-plot',[],{...PL},PC_EPICS);
    document.getElementById('physics-stats').textContent='';
}

function resizePhysics(){
    try{Plotly.Plots.resize('physics-plot');}catch(e){}
    try{Plotly.Plots.resize('moller-xy-plot');}catch(e){}
    try{Plotly.Plots.resize('hycal-xy-plot');}catch(e){}
}

function initPhysics(data){
    document.getElementById('physics-logz').onchange=()=>{plotEnergyAngle();plotMollerXY();plotHycalXY();};
    document.getElementById('physics-elastic').onchange=plotEnergyAngle;
}

// Theme flip — Plotly bakes THEME-derived colors (title/legend/colorbar
// titlefont, "ep elastic" line, etc.) into traces and layout at draw time,
// so a chrome-only relayout would leave stale colors behind. Replay the
// plot functions from cached *Data so titles/traces pick up live THEME.
if (typeof onThemeChange === 'function') {
    onThemeChange(() => {
        plotEnergyAngle();
        plotMollerXY();
        plotHycalXY();
    });
}
