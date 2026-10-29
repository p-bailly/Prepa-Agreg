// Accueil : série de jours, niveau atteint, compte à rebours de l'écrit, puis la liste de tous les chapitres, rangés par
// niveau (de la Terminale à l'agrégation) puis par thème, avec leur avancement. Toucher un chapitre rédigé l'ouvre.
import { useState } from 'react'
import { Volet } from '../composants/Blocs'
import { Cadenas, Pastille } from '../composants/Icones'
import { NIVEAUX, NOMS_NIVEAUX, type Niveau } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { etatAffiche } from '../repetition/etat'
import { joursAvantEcrit, niveauAtteint, progresNiveaux, serieJours } from '../stats/tableau'
import { etatChapitre, type EtatChapitre } from './Parcours'

/** « mar. 27 oct. » pour un jour civil `AAAA-MM-JJ` (midi UTC : le jour ne bascule dans aucun fuseau). */
const dateCourte = (jour: string) =>
  new Date(`${jour}T12:00:00Z`).toLocaleDateString('fr-FR', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' })

const PASTILLES = { acquis: 'plein', 'en-cours': 'moitie', 'a-venir': 'vide', vide: 'vide', masque: 'tiret' } as const
const NOMS_ETATS: Record<EtatChapitre, string> = {
  acquis: 'acquis',
  'en-cours': 'en cours',
  'a-venir': 'à découvrir',
  bloque: 'prérequis à voir d’abord',
  vide: 'pas encore rédigé',
  masque: 'masqué (sous le niveau de départ)',
}

export function Accueil({ surChapitre }: { surChapitre: (id: string) => void }) {
  const { contenu, progression, reglages, cartesParId, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const serie = serieJours(progression, jour)
  const niveau = niveauAtteint(progresNiveaux(contenu, progression, reglages))
  const nomNiveau = niveau ? NOMS_NIVEAUX[niveau] : '—'
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  // Ouverts d'emblée : les niveaux qui ont déjà des cartes rédigées.
  const [ouverts, setOuverts] = useState<Set<Niveau>>(() => new Set(contenu.cartes.map((k) => k.niveau)))
  const basculer = (n: Niveau, o: boolean) =>
    setOuverts((s) => {
      const t = new Set(s)
      if (o) t.add(n)
      else t.delete(n)
      return t
    })

  return (
    <main className="ecran">
      <div className="entete-accueil">
        <h1 className="titre-ecran">Prépa Agreg</h1>
        <span className="petit doux">{dateCourte(jour)}</span>
      </div>

      <div className="tuiles">
        <div className="tuile sable">
          <span className="nom">Série</span>
          <span className="valeur">{serie} j</span>
        </div>
        <div className="tuile menthe">
          <span className="nom">Niveau atteint</span>
          <span className={nomNiveau.length > 3 ? 'valeur longue' : 'valeur'}>{nomNiveau}</span>
        </div>
        {rebours !== null && (
          <div className="tuile rose">
            <span className="nom">Écrit</span>
            <span className={rebours > 0 ? 'valeur' : 'valeur longue'}>
              {rebours > 0 ? `${rebours} j` : rebours === 0 ? 'aujourd’hui' : 'passé'}
            </span>
          </div>
        )}
      </div>

      <h2 className="sous-titre">Chapitres</h2>
      {NIVEAUX.map((n) => {
        const themes = contenu.themes
          .map((t) => ({ theme: t, chapitres: t.chapitres.filter((ch) => ch.niveau === n) }))
          .filter((g) => g.chapitres.length)
        const cartes = contenu.cartes.filter((k) => k.niveau === n)
        const acquises = cartes.filter((k) => etatAffiche(k, progression, reglages) === 'acquise').length
        const nbChapitres = themes.reduce((s, g) => s + g.chapitres.length, 0)
        const marque = cartes.length ? `${acquises}/${cartes.length} acquises` : `${nbChapitres} chapitres à venir`
        return (
          <Volet key={n} titre={NOMS_NIVEAUX[n]} marque={marque} ouvert={ouverts.has(n)} surBascule={(o) => basculer(n, o)}>
            {themes.map(({ theme, chapitres }) => (
              <section key={theme.id} className="groupe-chapitres">
                <span className="etiquette">{theme.titre}</span>
                <ul className="liste-chapitres">
                  {chapitres.map((ch) => {
                    const e = etatChapitre(ch, contenu, progression, reglages, cartesParId)
                    const n = e.cartes.length
                    const inactif = e.etat === 'vide' || e.etat === 'masque'
                    return (
                      <li key={ch.id}>
                        <button className="chapitre-ligne" disabled={inactif} onClick={() => surChapitre(ch.id)}>
                          <span className="pastille" aria-hidden="true">
                            {e.etat === 'bloque' ? <Cadenas /> : <Pastille etat={PASTILLES[e.etat]} />}
                          </span>
                          <span>
                            {ch.titre}
                            <small>
                              {inactif ? `${NOMS_ETATS[e.etat]} · ${ch.prevues} notions prévues` : `${n} carte${n > 1 ? 's' : ''} · ${NOMS_ETATS[e.etat]}`}
                            </small>
                          </span>
                        </button>
                      </li>
                    )
                  })}
                </ul>
              </section>
            ))}
          </Volet>
        )
      })}
    </main>
  )
}
