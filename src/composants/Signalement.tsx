// « Signaler une erreur » : bouton discret et fenêtre de commentaire libre. Les signalements sont gardés dans la
// progression et inclus dans l'export (Réglages), pour être recollés à qui corrige le contenu.
import { useState } from 'react'
import { useDonnees } from '../etat/Donnees'
import { ajouterSignalement } from '../lib/progression'

export function Signaler({ cible, genre, ou }: { cible: string; genre: 'carte' | 'exercice'; ou?: string }) {
  const { changerProgression, horloge } = useDonnees()
  const [ouvert, setOuvert] = useState(false)
  const [texte, setTexte] = useState('')
  const [envoye, setEnvoye] = useState(false)
  const valider = () => {
    changerProgression((p) => ajouterSignalement(p, { date: horloge(), cible, genre, ou, commentaire: texte.trim() }))
    setOuvert(false)
    setTexte('')
    setEnvoye(true)
  }
  return (
    <>
      <button className="lien" onClick={() => setOuvert(true)}>
        {envoye ? 'Erreur signalée — en signaler une autre' : 'Signaler une erreur'}
      </button>
      {ouvert && (
        <div className="voile" onClick={() => setOuvert(false)}>
          <div className="fenetre" role="dialog" aria-label="Signaler une erreur" onClick={(e) => e.stopPropagation()}>
            <span className="etiquette">Signaler une erreur</span>
            <p className="petit doux">
              {genre === 'carte' ? 'Carte' : 'Exercice'} <code>{cible}</code>
              {ou ? ` · ${ou}` : ''}. Le signalement est gardé sur le téléphone et inclus dans l'export de la progression.
            </p>
            <textarea
              value={texte}
              onChange={(e) => setTexte(e.target.value)}
              placeholder="Ce qui semble faux, et pourquoi"
              autoFocus
            />
            <div className="deux-colonnes" style={{ marginTop: 12 }}>
              <button className="bouton" onClick={() => setOuvert(false)}>
                Annuler
              </button>
              <button className="bouton" disabled={!texte.trim()} onClick={valider}>
                Enregistrer
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  )
}
