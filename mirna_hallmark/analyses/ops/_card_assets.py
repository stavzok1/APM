"""CSS + JS for the card viewer — the browser half of `gen_cards.py`.

Split out for one reason: these are ~500 lines of a DIFFERENT language, and inlining them in the
generator made the Python unreadable. Nothing here is importable behaviour; it is two string
constants that `gen_cards.build()` interpolates into the page.

⭐ The design tokens are copied VERBATIM from `docs/derived/architecture.html` (`--ground --surface
--ink --muted --line --accent`, `prefers-color-scheme` + the `[data-theme]` overrides). Same
subproject, same look — a second palette would read as a second tool.

⛔ ZERO EXTERNAL REFERENCES. No CDN, no font, no image, no fetch. The pages must work from a USB
stick over `file://`, because that is how they reach a collaborator. `gen_cards --verify` greps the
output for `http://`, `https://` and `<script src` and fails on a hit.
"""
from __future__ import annotations

# --------------------------------------------------------------------------- #
CSS = """
:root{
  --ground:#f7f9fa; --surface:#ffffff; --surface-2:#eef2f4; --ink:#0f1720; --muted:#586a76;
  --line:#dde5e9; --accent:#0b7a8a; --accent-soft:#0b7a8a1a; --accent-ink:#075863;
  --s-strong:#1f9d6b; --s-mixed:#bd7d16; --s-open:#3f74d8; --s-design:#b0546a; --s-data:#7d6bc0;
  --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
}
@media (prefers-color-scheme:dark){:root{
  --ground:#0c1116; --surface:#131c23; --surface-2:#1b2730; --ink:#dbe6ed; --muted:#8496a3;
  --line:#25333d; --accent:#37c2d4; --accent-soft:#37c2d420; --accent-ink:#7fe0ec;
  --s-strong:#33c187; --s-mixed:#d99b2e; --s-open:#6f9bf0; --s-design:#d3778a; --s-data:#a898df;
}}
:root[data-theme="light"]{--ground:#f7f9fa;--surface:#ffffff;--surface-2:#eef2f4;--ink:#0f1720;--muted:#586a76;--line:#dde5e9;--accent:#0b7a8a;--accent-soft:#0b7a8a1a;--accent-ink:#075863;}
:root[data-theme="dark"]{--ground:#0c1116;--surface:#131c23;--surface-2:#1b2730;--ink:#dbe6ed;--muted:#8496a3;--line:#25333d;--accent:#37c2d4;--accent-soft:#37c2d420;--accent-ink:#7fe0ec;}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);line-height:1.5;-webkit-font-smoothing:antialiased}
a{color:var(--accent);text-decoration:none}
a:hover{text-decoration:underline}
.wrap{max-width:1240px;margin:0 auto;padding:clamp(18px,3.4vw,36px) clamp(14px,3vw,30px)}
header h1{font-size:clamp(1.35rem,3vw,2rem);font-weight:680;letter-spacing:-.02em;margin:0 0 .25em}
header .sub{color:var(--muted);max-width:74ch;margin:0}
.crumb{font-size:.85rem;margin-bottom:10px;color:var(--muted)}
.mono{font-family:var(--mono)}
.stats{display:flex;flex-wrap:wrap;gap:clamp(12px,3vw,30px);margin:18px 0 6px;padding:14px 0;border-top:1px solid var(--line);border-bottom:1px solid var(--line)}
.stat b{display:block;font-size:1.28rem;font-weight:660;letter-spacing:-.01em}
.stat span{color:var(--muted);font-size:.8rem}
.chip{display:inline-block;font-family:var(--mono);font-size:.72rem;padding:1px 6px;border-radius:5px;
  background:var(--surface-2);color:var(--muted);border:1px solid var(--line);white-space:nowrap}
.chip.rung{background:var(--accent-soft);color:var(--accent-ink);border-color:transparent}
.chip.agg{background:transparent;border-style:dashed}
.chip.warnchip{color:var(--s-mixed);border-color:var(--s-mixed);background:transparent;font-style:italic}
.searchbar{position:sticky;top:0;z-index:20;background:var(--ground);padding:12px 0 10px;
  display:flex;gap:10px;align-items:center;border-bottom:1px solid var(--line)}
input#q{flex:1;font:inherit;font-size:1rem;padding:9px 13px;border-radius:9px;border:1px solid var(--line);
  background:var(--surface);color:var(--ink)}
input#q:focus{outline:2px solid var(--accent-soft);border-color:var(--accent)}
#qn{color:var(--muted);font-size:.82rem;font-family:var(--mono);white-space:nowrap}
.grid{display:grid;grid-template-columns:minmax(230px,300px) 1fr;gap:22px;margin-top:18px;align-items:start}
@media(max-width:760px){.grid{grid-template-columns:1fr}}
.rail{position:sticky;top:64px;align-self:start;max-height:calc(100vh - 84px);overflow:auto;
  display:flex;flex-direction:column;gap:2px}
.hit{padding:7px 10px;border-radius:8px;cursor:pointer;border:1px solid transparent;font-size:.9rem;
  background:none;color:inherit;font-family:var(--mono);text-align:left;width:100%;
  overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.hit:hover{background:var(--surface-2)}
.hit.on{background:var(--accent-soft);color:var(--accent-ink);border-color:var(--accent)}
.railhead{font-size:.72rem;text-transform:uppercase;letter-spacing:.07em;color:var(--muted);
  padding:12px 10px 4px;font-weight:640}
main{min-width:0}
.card{background:var(--surface);border:1px solid var(--line);border-radius:13px;padding:clamp(14px,2.4vw,24px)}
.card h2{margin:0 0 2px;font-size:1.5rem;font-family:var(--mono);letter-spacing:-.01em;word-break:break-word}
.card .keyline{color:var(--muted);font-size:.85rem;margin:0 0 14px}
.blk{margin-top:20px;border-top:1px solid var(--line);padding-top:12px}
.blk>summary{cursor:pointer;font-weight:640;font-size:.95rem;list-style:none;display:flex;
  align-items:baseline;gap:9px;padding:2px 0}
.blk>summary::-webkit-details-marker{display:none}
.blk>summary::before{content:"\\25BE";color:var(--muted);font-size:.8em;transition:transform .15s}
.blk:not([open])>summary::before{transform:rotate(-90deg)}
.blk .n{color:var(--muted);font-weight:400;font-size:.8rem;font-family:var(--mono)}
.dom{font-size:.8rem;color:var(--muted);margin:7px 0 3px;padding:6px 10px;border-left:3px solid var(--line);
  background:var(--surface-2);border-radius:0 7px 7px 0}
.dom.warn{border-left-color:var(--s-mixed)}
.dom.stop{border-left-color:var(--s-design)}
table.vals{width:100%;border-collapse:collapse;font-size:.88rem}
table.vals td{padding:3px 8px 3px 0;vertical-align:baseline;border-bottom:1px solid var(--line)}
td.k{font-family:var(--mono);color:var(--muted);width:40%;word-break:break-all}
td.v{font-family:var(--mono);font-weight:600;white-space:nowrap}
td.v.na{font-weight:400;color:var(--muted);font-style:italic}
td.x{width:1%;white-space:nowrap;text-align:right}
svg.spark{vertical-align:middle;opacity:.62}
svg.spark rect{fill:var(--muted)}
svg.spark line{stroke:var(--accent);stroke-width:1.6}
.pctl{font-family:var(--mono);font-size:.72rem;color:var(--muted);margin-left:6px}
.charts{display:grid;grid-template-columns:repeat(auto-fit,minmax(215px,1fr));gap:16px;margin:16px 0 4px}
.ch{border:1px solid var(--line);border-radius:10px;padding:11px 13px;background:var(--surface-2)}
.ch h4{margin:0 0 3px;font-size:.82rem;font-weight:640}
.ch .cap{font-size:.72rem;color:var(--muted);margin:0 0 9px;line-height:1.35}
.bar{display:grid;grid-template-columns:74px 1fr 46px;gap:7px;align-items:center;font-size:.76rem;
  font-family:var(--mono);margin:3px 0}
.bar .t{height:9px;background:var(--line);border-radius:3px;overflow:hidden}
.bar .t i{display:block;height:100%;background:var(--accent);border-radius:3px}
.bar .val{text-align:right;color:var(--muted)}
.covwrap{display:flex;flex-wrap:wrap;gap:4px}
.cf{font-family:var(--mono);font-size:.68rem;padding:1px 5px;border-radius:4px;border:1px solid var(--line)}
.cf.on{background:var(--accent-soft);color:var(--accent-ink);border-color:transparent}
.cf.off{color:var(--muted);opacity:.62}
.dangle{color:var(--muted);font-style:italic}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(215px,1fr));gap:14px;margin:22px 0}
.tile{display:block;background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:16px;color:inherit}
.tile:hover{border-color:var(--accent);text-decoration:none}
.tile b{display:block;font-size:1.05rem;font-family:var(--mono);margin-bottom:3px}
.tile span{color:var(--muted);font-size:.8rem}
.note{font-size:.83rem;color:var(--muted);background:var(--surface-2);border:1px solid var(--line);
  border-radius:10px;padding:11px 14px;margin:16px 0;line-height:1.5}
.note.amber{border-color:var(--s-mixed);background:transparent}
footer{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);color:var(--muted);
  font-size:.74rem;font-family:var(--mono);word-break:break-all;line-height:1.7}
button.theme{position:fixed;right:14px;bottom:14px;z-index:40;border:1px solid var(--line);
  background:var(--surface);color:var(--muted);border-radius:9px;padding:6px 11px;cursor:pointer;font:inherit;font-size:.78rem}
.empty{color:var(--muted);padding:34px 6px;text-align:center;font-size:.92rem}
"""

# --------------------------------------------------------------------------- #
JS = r"""
const SEP='\u001f';   /* MUST equal gen_cards.SEP (unit separator) */
const P = JSON.parse(document.getElementById('payload').textContent);
const CACHE = {};
/* ⚠ Booleans survive the TSV -> pandas -> dictionary-encoder path as the STRINGS "True"/"False",
   NOT as JS booleans — verified against the encoder, not assumed. A `=== true` test silently never
   fires, which would have rendered every coverage flag as OFF and disabled the "not scanned" state
   entirely: the viewer would have looked perfectly fine and been wrong. Test through these two
   helpers, never against a literal. Neither coerces numbers — a measured 0 is not False. */
const isTrue  = v => v===true  || v==='True';
const isFalse = v => v===false || v==='False';

/* ---- decoder: the exact mirror of gen_cards.decode_column ------------------
   ⛔ The NaN mask is sacred here too. `x` lists rows that are blank DESPITE their
   group carrying a value; skipping that restore would show one arm's number on
   another arm's row.                                                          */
function dec(name){
  if(CACHE[name]) return CACHE[name];
  const e=P.d[name], t=e.t; let out;
  if(t==='n'){ out=e.v.split(',').map(x=> x===''?null:+x); }
  else if(t==='d'){ const L=e.l; out=e.i.split(',').map(x=> x===''?null:L[+x]); }
  else if(t==='s'){ out=e.v.split(SEP).map(x=> x===''?null:x); }
  else {
    const order=e.g.split(','), raw=(t==='gn'? e.v.split(',') : e.v.split(SEP));
    const gv=raw.map(x=> x===''?null:(t==='gn'? +x : x));
    const gi={}; order.forEach((g,i)=>gi[g]=i);
    out = dec(P.grouped[name]).map(g=>{ const j=gi[String(g)]; return j===undefined?null:gv[j]; });
    if(e.x) e.x.split(',').forEach(i=>{ if(i!=='') out[+i]=null; });
  }
  CACHE[name]=out; return out;
}
const esc = s => String(s).replace(/[&<>"']/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt = v => typeof v==='number' ? (Number.isInteger(v)? v.toLocaleString() :
                  (Math.abs(v)>=1e4||(Math.abs(v)<1e-3&&v!==0)? v.toExponential(2) : String(v))) : esc(v);

/* ---- search labels: one string per row, built from the key column(s) ------ */
const KEYS = P.key;
const LBL = (()=>{ const cols=KEYS.map(dec); const n=P.n, out=new Array(n);
  for(let i=0;i<n;i++) out[i]=cols.map(c=>c[i]===null?'':c[i]).join(' · ');
  return out; })();
const TOK = LBL.map(s=> s.toLowerCase().split(/[^a-z0-9.]+/i).filter(Boolean));

/* ⛔ NO fuzzy / edit-distance / subsequence matching. Gene symbols are dense
   (MIR21 vs MIR21HG; TP53 vs TP53BP1) and a confidently-WRONG entity rendered as
   a beautiful card is strictly worse than an empty result list. Five exact tiers.

   ⭐ Tier 3 (NAME BOUNDARY) exists because ties are broken by LENGTH, and without it
   `miR-21` returned hsa-miR-217, hsa-miR-2110, hsa-miR-2113 ABOVE hsa-miR-21-3p —
   measured, not hypothetical. The tokenizer splits on `-`, so a hyphenated query can
   never reach the token tiers and falls to plain substring, where "shortest first"
   ranks a DIFFERENT miRNA above the one asked for. Requiring the match to end on a
   non-alphanumeric separates `miR-21`-in-`miR-21-3p` from `miR-21`-in-`miR-217`. */
function btier(l,q){
  let i=l.indexOf(q), any=-1;
  while(i>=0){
    any=4;
    const a = i===0 || !/[a-z0-9]/.test(l[i-1]);
    const b = i+q.length===l.length || !/[a-z0-9]/.test(l[i+q.length]);
    if(a&&b) return 3;
    i=l.indexOf(q,i+1);
  }
  return any;
}
function rank(q){
  q=q.trim().toLowerCase(); if(!q) return [];
  const out=[];
  for(let i=0;i<LBL.length;i++){
    const l=LBL[i].toLowerCase(); let t;
    if(l===q) t=0;
    else if(TOK[i].indexOf(q)>=0) t=1;
    else if(TOK[i].some(x=>x.indexOf(q)===0)) t=2;
    else t=btier(l,q);
    if(t>=0){ out.push([t,LBL[i].length,i]); if(out.length>4000) break; }
  }
  out.sort((a,b)=> a[0]-b[0] || a[1]-b[1] || (LBL[a[2]]<LBL[b[2]]?-1:1));
  return out;
}
function colHits(q){
  q=q.trim().toLowerCase(); if(q.length<2) return [];
  return P.cols.filter(c=>c.toLowerCase().indexOf(q)>=0).slice(0,40);
}

/* ---- distribution strip --------------------------------------------------- */
function pctile(d,v){
  const q=d.q; if(v<=q[0]) return 0; if(v>=q[20]) return 100;
  for(let i=1;i<21;i++) if(v<=q[i]){ const lo=q[i-1],hi=q[i]; return Math.round((i-1+(hi>lo?(v-lo)/(hi-lo):0))*5); }
  return 100;
}
function spark(d,v){
  if(!d) return '';
  const h=d.h, mx=Math.max.apply(null,h), W=64, H=15, bw=W/h.length; let s='';
  for(let i=0;i<h.length;i++){ const bh=mx? Math.max(.8,H*h[i]/mx):.8;
    s+='<rect x="'+(i*bw).toFixed(2)+'" y="'+(H-bh).toFixed(2)+'" width="'+(bw-.55).toFixed(2)+'" height="'+bh.toFixed(2)+'"/>'; }
  if(v!=null && d.hi>d.lo){ const x=Math.min(W,Math.max(0,W*(v-d.lo)/(d.hi-d.lo)));
    s+='<line x1="'+x.toFixed(2)+'" y1="0" x2="'+x.toFixed(2)+'" y2="'+H+'"/>'; }
  return '<svg class="spark" viewBox="0 0 '+W+' '+H+'" width="'+W+'" height="'+H+'">'+s+'</svg>';
}

/* ---- the detail card ------------------------------------------------------ */
function render(i){
  if(i==null){ document.getElementById('detail').innerHTML =
    '<div class="empty">Search above, then pick a result.<br><span class="mono" style="font-size:.8em">'+
    P.n.toLocaleString()+' '+esc(P.card)+' rows · '+P.cols.length+' columns</span></div>'; return; }
  const row = {}; P.cols.forEach(c=> row[c]=dec(c)[i]);
  const title = KEYS.map(k=>row[k]).join(' · ');
  let h = '<div class="card"><h2>'+esc(title)+'</h2><p class="keyline">'+
     esc(P.card)+' card · key ('+KEYS.map(esc).join(', ')+') · row '+(i+1).toLocaleString()+
     ' of '+P.n.toLocaleString()+' · card rung <span class="chip rung">'+esc(P.own_rung)+'</span></p>';
  h += charts(row);
  for(const b of P.blocks){
    const cov = P.cov[b.k];
    /* ⭐ HIDE EMPTY BLOCKS. The median row fills a THIRD of its columns; without this the page is
       mostly blanks and the signal is unfindable. The count stays visible in the summary.
       ⛔ BUT NEVER HIDE AN *UNSCANNED* BLOCK. Hiding and unscanned are opposite messages, and the
       first version conflated them: 1,860 of 2,450 arms have no realization scan, and their
       Realization block VANISHED — the reader could not tell "this card has no such block" from
       "this scan never covered this arm", which is precisely the distinction the cov_ flags exist
       to carry. An unscanned block renders COLLAPSED and LABELLED instead: one line, still honest. */
    const filled = b.cols.filter(c=> row[c]!==null && row[c]!==undefined);
    const unscanned = !!(cov && isFalse(row[cov.flag]));
    if(!filled.length && !unscanned && b.k!=='coverage') continue;
    let warn='';
    const doms = b.cols.map(c=>(P.meta[c]||{}).dom||'').filter(Boolean);
    if(doms.some(d=>d.indexOf('⛔')>=0)) warn=' stop';
    else if(doms.some(d=>d.indexOf('⚠')>=0)) warn=' warn';
    h += '<details class="blk"'+((b.collapsed||!filled.length)?'':' open')+'><summary>'+esc(b.label)+
         ' <span class="n">'+filled.length+'/'+b.cols.length+'</span>'+
         (unscanned&&!filled.length ? ' <span class="chip warnchip">— not scanned ('+
            esc(cov.flag)+' is false)</span>' : '')+'</summary>';
    let last=null, rows='';
    for(const c of b.cols){
      const m=P.meta[c]||{}, dom=m.dom||'';
      if(dom && dom!==last){
        if(rows){ h+='<table class="vals">'+rows+'</table>'; rows=''; }
        const cls = dom.indexOf('⛔')>=0?' stop' : dom.indexOf('⚠')>=0?' warn':'';
        h += '<div class="dom'+cls+'">ⓘ '+esc(dom)+'</div>'; last=dom;
      }
      rows += valRow(c,row[c],m,cov,row);
    }
    if(rows) h+='<table class="vals">'+rows+'</table>';
    h += '</details>';
  }
  document.getElementById('detail').innerHTML = h+'</div>';
}

function valRow(c,v,m,cov,row){
  let cls='v', txt;
  if(v===null||v===undefined){
    cls='v na';
    /* ⭐ THREE DISTINCT EMPTY STATES, NEVER CONFLATED. A measured 0 is a finding; an UNSCANNED cell
       is an absence of evidence; a scanned-but-empty cell is neither. Rendering all three as a blank
       is the single most common way a reader invents a result. */
    txt = (cov && isFalse(row[cov.flag])) ? '— not scanned' : '— no value';
  } else txt = fmt(v);
  let chips='';
  /* the rung chip fires ONLY when the value is not at the card's own grain — otherwise 700 chips
     say the same thing and none of them is read */
  if(m.rung && m.rung!==P.own_rung && m.rung!=='key') chips+=' <span class="chip rung">⟨'+esc(m.rung)+'⟩</span>';
  if(m.agg) chips+=' <span class="chip agg">Σ '+esc(m.agg)+'</span>';
  let strip='';
  const d=P.dist[c];
  if(d && typeof v==='number'){ strip = spark(d,v)+'<span class="pctl">p'+pctile(d,v)+' of '+d.n.toLocaleString()+'</span>'; }
  /* cross-card link — resolved at BUILD time; dangling renders as text, never a link to a blank page */
  const L=P.links[c];
  if(L && v!==null && v!==undefined){
    const s=String(v);
    if(L.miss.indexOf(s)<0) txt='<a href="'+esc(L.to)+'.html#'+encodeURIComponent(s)+'">'+esc(s)+'</a>';
    else txt='<span class="dangle">'+esc(s)+' — '+(L.variant.indexOf(s)>=0?
        'on the '+esc(L.to)+' card under a different label form':'not on the '+esc(L.to)+' card')+'</span>';
  }
  return '<tr><td class="k">'+esc(c)+chips+'</td><td class="'+cls+'">'+txt+
         '</td><td class="x">'+strip+'</td></tr>';
}

/* ---- the three charts ----------------------------------------------------- */
function bars(items){       /* each row carries its OWN denominator; no shared axis */
  const mx=Math.max.apply(null,items.map(x=>x[1]).concat([1]));
  return items.map(x=>'<div class="bar"><span>'+esc(x[0])+'</span><span class="t"><i style="width:'+
    (100*x[1]/mx).toFixed(1)+'%"></i></span><span class="val">'+fmt(x[1])+'</span></div>').join('');
}
function charts(row){
  let out='';
  const ladder=['real_n_c05','real_n_c10','real_n_c15','real_n_c20','real_n_c30']
      .filter(c=>P.cols.indexOf(c)>=0 && row[c]!=null);
  if(ladder.length){
    out+='<div class="ch"><h4>Realization ladder</h4><p class="cap">Target genes whose coupling is at '+
      'least as strong as the cut. Nested by construction — the DECAY is the shape to read.</p>'+
      bars(ladder.map(c=>['≤ −'+(+c.slice(-2)/100).toFixed(2), row[c]]))+'</div>';
  }
  const cov=P.cols.filter(c=>c.indexOf('cov_')===0);
  if(cov.length){
    const on=cov.filter(c=>isTrue(row[c])).length;
    out+='<div class="ch"><h4>Coverage — '+on+'/'+cov.length+' scanned</h4>'+
      '<p class="cap">⚠ A sparse card is usually UNSCANNED, not uninteresting. These flags say which.</p>'+
      '<div class="covwrap">'+cov.map(c=>'<span class="cf '+(isTrue(row[c])?'on':'off')+'">'+
        esc(c.slice(4))+'</span>').join('')+'</div></div>';
  }
  const uni=[['scanMiR genome-wide','seq_n_genes_any'],['scanMiR site types','site_n_genes_canonical'],
             ['TargetScan','ts_n_genes'],['MANE 3′UTR seed','sdsz_N_7mer_plus']]
      .filter(u=>P.cols.indexOf(u[1])>=0 && row[u[1]]!=null);
  if(uni.length){
    /* ⛔ FOUR SEPARATE BARS, EACH WITH ITS OWN SCALE. One grouped chart with a shared axis is the
       VISUAL ACT of blending four incompatible targetome universes — different tools, different
       scopes, different denominators. The registry forbids it in words; this forbids it in pixels. */
    out+='<div class="ch"><h4>Targetome — four universes</h4><p class="cap">⛔ Never blended: '+
      'different tools, scopes and denominators. Each bar is scaled to ITSELF — lengths are '+
      'not comparable across rows.</p>'+
      uni.map(u=>{const d=P.dist[u[1]]; const v=row[u[1]];
        return '<div class="bar"><span>'+esc(u[0])+'</span><span class="t"><i style="width:'+
        (d&&d.hi>d.lo? (100*Math.min(1,(v-d.lo)/(d.hi-d.lo))).toFixed(1):'0')+
        '%"></i></span><span class="val">'+fmt(v)+'</span></div>';}).join('')+'</div>';
  }
  return out? '<div class="charts">'+out+'</div>' : '';
}

/* ---- wiring --------------------------------------------------------------- */
let hits=[], sel=-1;
function draw(){
  const q=document.getElementById('q').value;
  hits=rank(q); const cols=colHits(q);
  document.getElementById('qn').textContent = q.trim()? hits.length.toLocaleString()+' row'+(hits.length===1?'':'s') : P.n.toLocaleString()+' rows';
  let h = hits.slice(0,300).map((x,j)=>'<button class="hit'+(j===sel?' on':'')+'" data-i="'+x[2]+'" data-j="'+j+'">'+esc(LBL[x[2]])+'</button>').join('');
  if(hits.length>300) h+='<div class="railhead">+'+(hits.length-300).toLocaleString()+' more — refine</div>';
  if(cols.length) h+='<div class="railhead">columns on this card</div>'+
     cols.map(c=>'<button class="hit" data-col="'+esc(c)+'">'+esc(c)+'</button>').join('');
  document.getElementById('rail').innerHTML = h || '<div class="railhead">no match</div>';
}
document.getElementById('rail').addEventListener('click',e=>{
  const b=e.target.closest('.hit'); if(!b) return;
  if(b.dataset.col!==undefined){ document.getElementById('q').value=b.dataset.col; draw(); return; }
  sel=+b.dataset.j; location.hash=encodeURIComponent(LBL[+b.dataset.i]); draw(); render(+b.dataset.i);
});
document.getElementById('q').addEventListener('input',()=>{sel=-1;draw();});
document.addEventListener('keydown',e=>{
  const q=document.getElementById('q');
  if((e.key==='/'||(e.key==='k'&&(e.metaKey||e.ctrlKey))) && document.activeElement!==q){e.preventDefault();q.focus();q.select();return;}
  if(e.key==='Escape'){q.value='';sel=-1;draw();q.blur();return;}
  if(!hits.length) return;
  if(e.key==='ArrowDown'||e.key==='ArrowUp'){
    e.preventDefault(); sel=Math.max(0,Math.min(hits.length-1,sel+(e.key==='ArrowDown'?1:-1)));
    draw(); const el=document.querySelector('.hit.on'); if(el) el.scrollIntoView({block:'nearest'});
  } else if(e.key==='Enter'&&sel>=0){ const i=hits[sel][2]; location.hash=encodeURIComponent(LBL[i]); render(i); }
});
function fromHash(){
  const k=decodeURIComponent(location.hash.slice(1)); if(!k) { render(null); return; }
  let i=LBL.indexOf(k);
  if(i<0){ const lo=k.toLowerCase(); i=LBL.findIndex(s=>s.toLowerCase()===lo); }
  if(i<0){ document.getElementById('q').value=k; draw(); render(null); return; }
  document.getElementById('q').value=k; draw(); render(i);
}
window.addEventListener('hashchange',fromHash);
document.querySelector('button.theme').addEventListener('click',()=>{
  const r=document.documentElement;
  const cur=r.getAttribute('data-theme')||(matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light');
  r.setAttribute('data-theme',cur==='dark'?'light':'dark');
});
draw(); fromHash();
"""

# --------------------------------------------------------------------------- #
INDEX_JS = r"""
const P = JSON.parse(document.getElementById('payload').textContent);
const esc = s => String(s).replace(/[&<>"']/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const ROWS=[];   /* [label, card] over every entity on every card */
for(const c of P.order) for(const k of P.cards[c].keys) ROWS.push([k,c]);
const TOK = ROWS.map(r=> r[0].toLowerCase().split(/[^a-z0-9.]+/i).filter(Boolean));
/* same five-tier matcher as the card pages — see the note on btier() in JS. Duplicated rather than
   shared because each page is standalone; `--verify-js` asserts both give the same ranking. */
function btier(l,q){
  let i=l.indexOf(q), any=-1;
  while(i>=0){ any=4;
    const a = i===0 || !/[a-z0-9]/.test(l[i-1]);
    const b = i+q.length===l.length || !/[a-z0-9]/.test(l[i+q.length]);
    if(a&&b) return 3;
    i=l.indexOf(q,i+1); }
  return any;
}

function draw(){
  const q=document.getElementById('q').value.trim().toLowerCase();
  if(!q){ document.getElementById('rail').innerHTML=''; document.getElementById('qn').textContent=
      P.total.toLocaleString()+' rows · '+P.ncols+' columns'; return; }
  const out=[];
  for(let i=0;i<ROWS.length;i++){
    const l=ROWS[i][0].toLowerCase(); let t;
    if(l===q) t=0; else if(TOK[i].indexOf(q)>=0) t=1;
    else if(TOK[i].some(x=>x.indexOf(q)===0)) t=2; else t=btier(l,q);
    if(t>=0){ out.push([t,ROWS[i][0].length,i]); if(out.length>3000) break; }
  }
  out.sort((a,b)=>a[0]-b[0]||a[1]-b[1]);
  /* ⭐ "where does this column live" — the question a 699-column system makes unanswerable
     without a tool, and the reason column names are searchable at all. */
  const cols=[]; if(q.length>=2) for(const c of P.order)
    for(const n of P.cards[c].cols) if(n.toLowerCase().indexOf(q)>=0) cols.push([n,c]);
  document.getElementById('qn').textContent = out.length.toLocaleString()+' row'+(out.length===1?'':'s')+
      (cols.length? ' · '+cols.length+' column'+(cols.length===1?'':'s'):'');
  let h = out.slice(0,220).map(x=>'<a class="hit" href="'+ROWS[x[2]][1]+'.html#'+
      encodeURIComponent(ROWS[x[2]][0])+'">'+esc(ROWS[x[2]][0])+
      ' <span class="chip">'+ROWS[x[2]][1]+'</span></a>').join('');
  if(out.length>220) h+='<div class="railhead">+'+(out.length-220).toLocaleString()+' more — refine</div>';
  if(cols.length) h+='<div class="railhead">columns — which card carries them</div>'+
     cols.slice(0,60).map(x=>'<div class="hit">'+esc(x[0])+' <span class="chip">'+x[1]+'</span></div>').join('');
  document.getElementById('rail').innerHTML=h||'<div class="railhead">no match</div>';
}
document.getElementById('q').addEventListener('input',draw);
document.addEventListener('keydown',e=>{
  const q=document.getElementById('q');
  if((e.key==='/'||(e.key==='k'&&(e.metaKey||e.ctrlKey)))&&document.activeElement!==q){e.preventDefault();q.focus();q.select();}
  if(e.key==='Escape'){q.value='';draw();q.blur();}
});
document.querySelector('button.theme').addEventListener('click',()=>{
  const r=document.documentElement;
  const cur=r.getAttribute('data-theme')||(matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light');
  r.setAttribute('data-theme',cur==='dark'?'light':'dark');
});
draw();
"""
