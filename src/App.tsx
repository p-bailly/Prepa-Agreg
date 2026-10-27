import { useState } from 'react'

type Onglet = 'accueil' | 'parcours' | 'stats' | 'reglages'

const ONGLETS: { id: Onglet; nom: string }[] = [
  { id: 'accueil', nom: 'Accueil' },
  { id: 'parcours', nom: 'Parcours' },
  { id: 'stats', nom: 'Stats' },
  { id: 'reglages', nom: 'Réglages' },
]

export default function App() {
  const [onglet, setOnglet] = useState<Onglet>('accueil')
  const actuel = ONGLETS.find((o) => o.id === onglet)!
  return (
    <div className="appli avec-onglets">
      <main className="ecran">
        <h1 className="titre-ecran">{onglet === 'accueil' ? 'Prépa Agreg' : actuel.nom}</h1>
        <span className="etiquette">Bientôt</span>
        <div className="bloc">
          <p>Écran en construction.</p>
        </div>
      </main>
      <nav className="onglets" aria-label="Onglets">
        {ONGLETS.map((o) => (
          <button key={o.id} aria-current={o.id === onglet ? 'page' : undefined} onClick={() => setOnglet(o.id)}>
            {o.nom}
          </button>
        ))}
      </nav>
    </div>
  )
}
