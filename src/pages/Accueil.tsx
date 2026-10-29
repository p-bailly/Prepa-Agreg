// Accueil : compte à rebours de l'écrit (si une date est choisie), puis la liste de tous les chapitres, numérotés, un
// volet par niveau (de la Terminale à l'agrégation) et groupés par thème, avec leur avancement. Toucher un chapitre
// rédigé l'ouvre (notions, exercices et problèmes, QCM).
import { useMemo, useState } from 'react'
import { avancementChapitre, chapitresNumerotes, notionsDe } from '../chapitres/plan'
import { Volet } from '../composants/Blocs'
import { Ligne, NOMS_AVANCEMENT } from '../composants/Lignes'
import { NIVEAUX, NOMS_NIVEAUX, type Niveau } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { joursAvantEcrit } from '../lib/temps'

/** « mar. 27 oct. » pour un jour civil `AAAA-MM-JJ` (midi UTC : le jour ne bascule dans aucun fuseau). */
const dateCourte = (jour: string) =>
  new Date(`${jour}T12:00:00Z`).toLocaleDateString('fr-FR', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' })

export function Accueil({ surChapitre }: { surChapitre: (id: string) => void }) {
  const { contenu, progression, reglages, exercicesParCarte, aujourdhui } = useDonnees()
  const jour = aujourdhui()
  const rebours = joursAvantEcrit(reglages.dateEcrit, jour)
  const chapitres = useMemo(() => chapitresNumerotes(contenu), [contenu])
  // Ouverts d'emblée : les niveaux qui ont déjà des chapitres rédigés.
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

      {rebours !== null && (
        <div className="tuiles">
          <div className="tuile rose">
            <span className="nom">Écrit</span>
            <span className={rebours > 0 ? 'valeur' : 'valeur longue'}>
              {rebours > 0 ? `${rebours} j` : rebours === 0 ? 'aujourd’hui' : 'passé'}
            </span>
          </div>
        </div>
      )}

      {NIVEAUX.map((n) => {
        const duNiveau = chapitres.filter((x) => x.chapitre.niveau === n)
        const themes = [...new Set(duNiveau.map((x) => x.theme))]
        const rediges = duNiveau.filter((x) => notionsDe(contenu, x.chapitre.id).length).length
        const premier = duNiveau[0]?.numero
        const dernier = duNiveau[duNiveau.length - 1]?.numero
        const marque = rediges ? `ch. ${premier} à ${dernier} · ${rediges}/${duNiveau.length} rédigés` : `ch. ${premier} à ${dernier} · à venir`
        return (
          <Volet key={n} titre={NOMS_NIVEAUX[n]} marque={marque} ouvert={ouverts.has(n)} surBascule={(o) => basculer(n, o)}>
            {themes.map((theme) => (
              <section key={theme.id} className="groupe-chapitres">
                <span className="etiquette">{theme.titre}</span>
                <ul className="liste">
                  {duNiveau
                    .filter((x) => x.theme === theme)
                    .map(({ numero, chapitre }) => {
                      const a = avancementChapitre(contenu, chapitre.id, progression)
                      const notions = notionsDe(contenu, chapitre.id)
                      const nbEx = notions.reduce((s, k) => s + (exercicesParCarte.get(k.id)?.length ?? 0), 0)
                      return (
                        <li key={chapitre.id}>
                          <Ligne
                            avancement={a}
                            numero={String(numero)}
                            titre={chapitre.titre}
                            inactif={a === 'vide'}
                            dessous={
                              a === 'vide'
                                ? `${NOMS_AVANCEMENT.vide} · ${chapitre.prevues} notions prévues`
                                : `${notions.length} notions · ${nbEx} exercices · QCM · ${NOMS_AVANCEMENT[a]}`
                            }
                            surClic={() => surChapitre(chapitre.id)}
                          />
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
