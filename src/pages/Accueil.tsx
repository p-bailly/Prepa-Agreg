// Accueil : la liste de tous les chapitres, numérotés, un volet par niveau (de la Terminale à l'agrégation interne) et
// groupés par thème (icône et couleur du thème), en cartes sur une colonne au téléphone, sur deux ou trois sur grand
// écran. Toucher un chapitre rédigé l'ouvre (notions, exercices
// et problèmes, QCM). En tête de chaque niveau rédigé : l'avancement du niveau, un QCM et un vrai-faux sur tout le niveau.
import { useMemo, useState } from 'react'
import { avancementChapitre, chapitresNumerotes, meilleurQcm, notionsDe } from '../chapitres/plan'
import { cleQcm, NB_QUESTIONS, questionsNiveau } from '../chapitres/qcm'
import { Volet } from '../composants/Blocs'
import { IconeTheme } from '../composants/Icones'
import { Ligne, NOMS_AVANCEMENT } from '../composants/Lignes'
import { NIVEAUX, NOMS_NIVEAUX, type Niveau } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'

/** « mar. 27 oct. » pour un jour civil `AAAA-MM-JJ` (midi UTC : le jour ne bascule dans aucun fuseau). */
const dateCourte = (jour: string) =>
  new Date(`${jour}T12:00:00Z`).toLocaleDateString('fr-FR', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' })

export function Accueil({
  surChapitre,
  surQcmNiveau,
}: {
  surChapitre: (id: string) => void
  surQcmNiveau: (niveau: Niveau, sorte: 'qcm' | 'vrai-faux') => void
}) {
  const { contenu, progression, exercicesParCarte, aujourdhui } = useDonnees()
  const jour = aujourdhui()
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

      {NIVEAUX.map((n) => {
        const duNiveau = chapitres.filter((x) => x.chapitre.niveau === n)
        const themes = [...new Set(duNiveau.map((x) => x.theme))]
        const rediges = duNiveau.filter((x) => notionsDe(contenu, x.chapitre.id).length).length
        const faits = duNiveau.filter((x) => avancementChapitre(contenu, x.chapitre.id, progression) === 'fait').length
        const premier = duNiveau[0]?.numero
        const dernier = duNiveau[duNiveau.length - 1]?.numero
        const marque = rediges ? `ch. ${premier} à ${dernier} · ${rediges}/${duNiveau.length} rédigés` : `ch. ${premier} à ${dernier} · à venir`
        return (
          <Volet key={n} titre={NOMS_NIVEAUX[n]} marque={marque} ouvert={ouverts.has(n)} surBascule={(o) => basculer(n, o)}>
            {rediges > 0 && (
              <div className="barre-niveau">
                <span className="piste" aria-hidden="true">
                  <span style={{ width: `${(100 * faits) / duNiveau.length}%` }} />
                </span>
                {faits} {faits > 1 ? 'chapitres faits' : 'chapitre fait'} sur {duNiveau.length}
              </div>
            )}
            {rediges > 0 && (
              <section className="groupe-chapitres">
                <span className="etiquette">Tout le niveau</span>
                <div className="niveau-entier">
                  {(['qcm', 'vrai-faux'] as const).map((sorte) => {
                    const total = questionsNiveau(contenu, n, sorte).length
                    const meilleur = meilleurQcm(progression, cleQcm({ niveau: n, sorte }))
                    return (
                      <button key={sorte} className="bouton" disabled={!total} onClick={() => surQcmNiveau(n, sorte)}>
                        {sorte === 'qcm' ? 'QCM du niveau' : 'Vrai ou faux du niveau'}
                        <small>
                          {Math.min(NB_QUESTIONS.niveau, total)} questions parmi {total}
                          {meilleur ? ` · meilleur ${meilleur.bonnes}/${meilleur.total}` : ''}
                        </small>
                      </button>
                    )
                  })}
                </div>
              </section>
            )}
            <div className="themes-niveau">
              {themes.map((theme) => (
                <section key={theme.id} className={`groupe-chapitres th-${theme.id}`}>
                  <span className="etiquette-theme">
                    <span className="rond">
                      <IconeTheme theme={theme.id} />
                    </span>
                    {theme.titre}
                  </span>
                  <ul className="grille-chapitres">
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
            </div>
          </Volet>
        )
      })}
    </main>
  )
}
