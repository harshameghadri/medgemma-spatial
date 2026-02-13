import Head from 'next/head';
import '../styles/globals.css';

const HF_SPACE_URL = 'https://harshameghadri-medgemma-spatial.hf.space';
const GITHUB_URL = 'https://github.com/harshameghadri/medgemma-spatial';
const KAGGLE_URL = 'https://www.kaggle.com/code/harshameghadri/medgemma-spatial-transcriptomics-analysis';
const NBVIEWER_URL = 'https://nbviewer.org/github/harshameghadri/medgemma-spatial/blob/main/notebooks/05_whitepaper.ipynb';

const metrics = [
  { label: 'Tissue Spots Analysed', value: '~5K', unit: 'per sample', color: '#3b82f6' },
  { label: 'Cell Types Annotated', value: '13', unit: 'compartments', color: '#8b5cf6' },
  { label: 'Report Quality Metrics', value: '5', unit: 'benchmark axes', color: '#22c55e' },
];

const techStack = [
  'MedGemma 4B-it', 'QLoRA', 'Scanpy', 'Squidpy',
  'CellTypist', 'Streamlit', 'HuggingFace', 'Python 3.10',
];

const pipeline = [
  { step: '01', title: 'QC & Preprocessing', desc: 'Filter spots, normalise expression, select 2K HVGs, Leiden clustering' },
  { step: '02', title: 'Tissue Annotation', desc: 'Tissue-agnostic 2-tier: z-score marker panels + selective CellTypist immune pass' },
  { step: '03', title: 'Spatial Statistics', desc: "Moran's I (100 permutations), bootstrap entropy CIs, multi-scale neighbourhood enrichment" },
  { step: '04', title: 'Feature JSON', desc: 'Structured spatial features: annotation + heterogeneity + uncertainty' },
  { step: '05', title: 'MedGemma Report', desc: 'MoA-focused 5-section clinical report: TME → mechanism → research → therapy → caveats' },
];

export default function Home() {
  return (
    <>
      <Head>
        <title>MedGemma Spatial Transcriptomics | Portfolio</title>
        <meta name="description" content="Automated mechanism-of-action synthesis from spatial transcriptomics using MedGemma + QLoRA fine-tuning." />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <link rel="icon" href="/favicon.ico" />
      </Head>

      <div style={{ maxWidth: 1100, margin: '0 auto', padding: '0 1.5rem' }}>

        {/* Nav */}
        <nav style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', padding: '1.5rem 0', borderBottom: '1px solid #1e2d45' }}>
          <span style={{ fontWeight: 700, fontSize: '1.1rem', color: '#3b82f6' }}>MedGemma Spatial</span>
          <div style={{ display: 'flex', gap: '1.5rem', fontSize: '0.9rem' }}>
            <a href={GITHUB_URL} target="_blank" rel="noopener noreferrer">GitHub</a>
            <a href={KAGGLE_URL} target="_blank" rel="noopener noreferrer">Kaggle</a>
            <a href={NBVIEWER_URL} target="_blank" rel="noopener noreferrer">Whitepaper</a>
          </div>
        </nav>

        {/* Hero */}
        <section style={{ padding: '5rem 0 3rem', textAlign: 'center' }}>
          <div style={{ display: 'inline-block', background: 'rgba(59,130,246,0.1)', border: '1px solid rgba(59,130,246,0.3)', borderRadius: 20, padding: '0.35rem 1rem', fontSize: '0.8rem', color: '#93c5fd', marginBottom: '1.5rem', letterSpacing: '0.05em', textTransform: 'uppercase' }}>
            Google MedGemma Impact Challenge 2026
          </div>
          <h1 style={{ fontSize: 'clamp(2rem, 5vw, 3.5rem)', fontWeight: 800, lineHeight: 1.15, marginBottom: '1.5rem', background: 'linear-gradient(135deg, #e2e8f0 0%, #93c5fd 100%)', WebkitBackgroundClip: 'text', WebkitTextFillColor: 'transparent' }}>
            Automated MoA Synthesis<br />from Spatial Transcriptomics
          </h1>
          <p style={{ fontSize: '1.15rem', color: '#94a3b8', maxWidth: 640, margin: '0 auto 2.5rem', lineHeight: 1.7 }}>
            A full-stack AI pipeline that converts 10x Visium data into mechanism-of-action focused clinical pathology reports using MedGemma + QLoRA fine-tuning.
          </p>
          <div style={{ display: 'flex', gap: '1rem', justifyContent: 'center', flexWrap: 'wrap' }}>
            <a href="#demo" style={{ background: '#3b82f6', color: '#fff', padding: '0.75rem 2rem', borderRadius: 8, fontWeight: 600, textDecoration: 'none' }}>
              Live Demo ↓
            </a>
            <a href={NBVIEWER_URL} target="_blank" rel="noopener noreferrer" style={{ background: 'transparent', color: '#e2e8f0', padding: '0.75rem 2rem', borderRadius: 8, fontWeight: 600, border: '1px solid #1e2d45', textDecoration: 'none' }}>
              Read Whitepaper
            </a>
          </div>
        </section>

        {/* Metrics */}
        <section style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem', margin: '2rem 0' }}>
          {metrics.map(m => (
            <div key={m.label} style={{ background: '#111827', border: '1px solid #1e2d45', borderRadius: 12, padding: '1.5rem', textAlign: 'center' }}>
              <div style={{ fontSize: '2.5rem', fontWeight: 800, color: m.color }}>{m.value}</div>
              <div style={{ fontSize: '0.75rem', color: '#94a3b8', textTransform: 'uppercase', letterSpacing: '0.05em' }}>{m.unit}</div>
              <div style={{ fontSize: '0.9rem', color: '#cbd5e1', marginTop: '0.25rem' }}>{m.label}</div>
            </div>
          ))}
        </section>

        {/* Pipeline */}
        <section style={{ margin: '4rem 0' }}>
          <h2 style={{ fontSize: '1.75rem', fontWeight: 700, marginBottom: '2rem', textAlign: 'center' }}>5-Stage Pipeline</h2>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
            {pipeline.map((p, i) => (
              <div key={p.step} style={{ display: 'flex', gap: '1.25rem', background: '#111827', border: '1px solid #1e2d45', borderRadius: 10, padding: '1.25rem 1.5rem', alignItems: 'flex-start' }}>
                <span style={{ fontSize: '0.8rem', fontWeight: 700, color: '#3b82f6', minWidth: 28, marginTop: 3 }}>{p.step}</span>
                <div>
                  <div style={{ fontWeight: 600, marginBottom: '0.25rem' }}>{p.title}</div>
                  <div style={{ fontSize: '0.875rem', color: '#94a3b8' }}>{p.desc}</div>
                </div>
              </div>
            ))}
          </div>
        </section>

        {/* Tech Stack */}
        <section style={{ margin: '3rem 0', textAlign: 'center' }}>
          <h2 style={{ fontSize: '1.4rem', fontWeight: 600, marginBottom: '1.25rem', color: '#94a3b8' }}>Built With</h2>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.6rem', justifyContent: 'center' }}>
            {techStack.map(t => (
              <span key={t} style={{ background: '#1e2d45', border: '1px solid #2d3f5a', borderRadius: 6, padding: '0.35rem 0.9rem', fontSize: '0.85rem', color: '#93c5fd' }}>{t}</span>
            ))}
          </div>
        </section>

        {/* Live Demo */}
        <section id="demo" style={{ margin: '4rem 0' }}>
          <h2 style={{ fontSize: '1.75rem', fontWeight: 700, marginBottom: '0.75rem', textAlign: 'center' }}>Live Demo</h2>
          <p style={{ textAlign: 'center', color: '#94a3b8', marginBottom: '1.5rem', fontSize: '0.9rem' }}>
            Upload a 10x Visium .h5ad file to get an automated clinical pathology report.
            Hosted on HuggingFace Spaces.
          </p>
          <div style={{ border: '1px solid #1e2d45', borderRadius: 12, overflow: 'hidden', background: '#111827' }}>
            <div style={{ padding: '0.75rem 1rem', borderBottom: '1px solid #1e2d45', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
              <div style={{ width: 10, height: 10, borderRadius: '50%', background: '#ef4444' }} />
              <div style={{ width: 10, height: 10, borderRadius: '50%', background: '#f59e0b' }} />
              <div style={{ width: 10, height: 10, borderRadius: '50%', background: '#22c55e' }} />
              <span style={{ marginLeft: '0.5rem', fontSize: '0.8rem', color: '#64748b' }}>{HF_SPACE_URL}</span>
            </div>
            <iframe
              src={HF_SPACE_URL}
              width="100%"
              height="700"
              style={{ border: 'none', display: 'block' }}
              title="MedGemma Spatial Transcriptomics Demo"
              loading="lazy"
            />
          </div>
          <p style={{ textAlign: 'center', marginTop: '0.75rem', fontSize: '0.8rem', color: '#475569' }}>
            App may take 30-60s to wake up on first load.{' '}
            <a href={HF_SPACE_URL} target="_blank" rel="noopener noreferrer">Open in full screen →</a>
          </p>
        </section>

        {/* Links */}
        <section style={{ margin: '3rem 0', display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))', gap: '1rem' }}>
          {[
            { title: '📓 Kaggle Notebook', desc: 'Full analysis pipeline on public Visium dataset', href: KAGGLE_URL },
            { title: '📄 Whitepaper Notebook', desc: '2-page technical narrative with benchmarks', href: NBVIEWER_URL },
            { title: '💻 GitHub Repository', desc: 'Source code, scripts, training pipeline', href: GITHUB_URL },
          ].map(l => (
            <a key={l.title} href={l.href} target="_blank" rel="noopener noreferrer" style={{ background: '#111827', border: '1px solid #1e2d45', borderRadius: 10, padding: '1.5rem', display: 'block', textDecoration: 'none', transition: 'border-color 0.2s' }}>
              <div style={{ fontWeight: 600, marginBottom: '0.5rem', color: '#e2e8f0' }}>{l.title}</div>
              <div style={{ fontSize: '0.85rem', color: '#94a3b8' }}>{l.desc}</div>
            </a>
          ))}
        </section>

        {/* Footer */}
        <footer style={{ borderTop: '1px solid #1e2d45', padding: '2rem 0', textAlign: 'center', color: '#475569', fontSize: '0.85rem' }}>
          <p>Built by Sriharsha Meghadri · Google MedGemma Impact Challenge 2026</p>
          <p style={{ marginTop: '0.5rem' }}>
            <a href={GITHUB_URL} target="_blank" rel="noopener noreferrer">GitHub</a>
            {' · '}
            <a href={KAGGLE_URL} target="_blank" rel="noopener noreferrer">Kaggle</a>
            {' · '}
            <a href={HF_SPACE_URL} target="_blank" rel="noopener noreferrer">HuggingFace</a>
          </p>
        </footer>
      </div>
    </>
  );
}
