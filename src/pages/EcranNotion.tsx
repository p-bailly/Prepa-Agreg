// Notion d'un chapitre : deux onglets. Cours : texte de cours, énoncé encadré, exemple, démonstration, contre-exemples,
// pièges (vrai ou faux), et un bouton pour mémoriser l'énoncé avec la carte de révision. Exercices : rangés par niveau de
// difficulté (1 application directe … 5 difficile). Ouvrir le cours marque la notion comme apprise.
import { useEffect, useMemo, useState } from 'react'
import { chapitresNumerotes, notionsDe } from '../chapitres/plan'
import { BandeauDoute, Encadre, Segments, Volet } from '../composants/Blocs'
import { LigneExercice } from '../composants/Lignes'
import { Signaler } from '../composants/Signalement'
import { Texte } from '../composants/Texte'
import type { Carte } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { etatExercice } from '../exercices/niveaux'
import { apprendre } from '../lib/progression'

type Onglet = 'cours' | 'exercices'

/** Niveaux de difficulté, voir CLAUDE.md, « Exercice ». */
const NIVEAUX_DIFFICULTE = [
  'application directe',
  'calcul ou application standard',
  'classique en plusieurs étapes',
  'question d’écrit',
  'difficile',
]

export function EcranNotion({
  carte,
  surExercice,
  surMemoriser,
  surRetour,
}: {
  carte: Carte
  surExercice: (id: string) => void
  surMemoriser: () => void
  surRetour: () => void
}) {
  const { contenu, progression, changerProgression, exercicesParCarte, aujourdhui } = useDonnees()
  const [onglet, setOnglet] = useState<Onglet>('cours')
  const entree = useMemo(() => chapitresNumerotes(contenu).find((x) => x.chapitre.id === carte.chapitre), [contenu, carte.chapitre])
  const rang = useMemo(() => notionsDe(contenu, carte.chapitre).findIndex((k) => k.id === carte.id) + 1, [contenu, carte])
  const exercices = exercicesParCarte.get(carte.id) ?? []

  useEffect(() => {
    if (onglet === 'cours') changerProgression((p) => apprendre(p, [carte.id], aujourdhui()))
  }, [onglet, carte.id, changerProgression, aujourdhui])

  return (
    <main className={`ecran lecture th-${carte.theme}`}>
      <button className="retour" onClick={surRetour}>
        ← Chapitre {entree?.numero}
      </button>
      <span className="sur-titre-ecran">
        Notion {entree?.numero}.{rang} · {entree?.chapitre.titre}
      </span>
      <h1 className="titre-ecran">
        <Texte texte={carte.titre} enLigne />
      </h1>
      <Segments
        options={[
          { id: 'cours', nom: 'Cours' },
          { id: 'exercices', nom: 'Exercices' },
        ]}
        valeur={onglet}
        surChange={setOnglet}
      />

      {onglet === 'cours' && (
        <>
          <BandeauDoute doute={carte.doute} />
          {carte.cours && <Texte texte={carte.cours} />}
          <Encadre genre={carte.type}>
            <Texte texte={carte.enonce} />
          </Encadre>
          {carte.exemple && (
            <Encadre genre="exemple">
              <Texte texte={carte.exemple} />
            </Encadre>
          )}
          {carte.demo && (
            <Volet titre="Démonstration">
              <div style={{ paddingTop: 12 }}>
                <Texte texte={carte.demo} />
              </div>
            </Volet>
          )}
          {carte.contre_exemple && (
            <Volet titre="Pourquoi ces hypothèses">
              <Encadre genre="remarque" nom="Contre-exemples">
                <Texte texte={carte.contre_exemple} />
              </Encadre>
            </Volet>
          )}
          {!!carte.affirmations?.length && (
            <Volet titre="Pièges à éviter">
              <div style={{ paddingTop: 12 }}>
                {carte.affirmations.map((a, i) => (
                  <Encadre key={i} genre={a.vrai ? 'definition' : 'theoreme'} nom={a.vrai ? 'Vrai' : 'Faux'}>
                    <Texte texte={a.texte} />
                    {a.vrai && a.justification && <Texte texte={a.justification} />}
                    {!a.vrai && a.contre_exemple && (
                      <>
                        <p className="petit doux" style={{ marginTop: 8 }}>
                          Contre-exemple :
                        </p>
                        <Texte texte={a.contre_exemple} />
                      </>
                    )}
                  </Encadre>
                ))}
              </div>
            </Volet>
          )}
          <button className="bouton" onClick={surMemoriser}>
            Mémoriser l’énoncé (carte de révision)
          </button>
          <button className="bouton" onClick={() => setOnglet('exercices')}>
            Passer aux exercices
          </button>
          <Signaler cible={carte.id} genre="carte" ou="cours" />
        </>
      )}

      {onglet === 'exercices' &&
        (exercices.length ? (
          [1, 2, 3, 4, 5].map((d) => {
            const l = exercices.filter((e) => e.difficulte === d)
            if (!l.length) return null
            return (
              <section key={d}>
                <h2 className="sous-titre">
                  Niveau {d} <span className="doux petit">· {NIVEAUX_DIFFICULTE[d - 1]}</span>
                </h2>
                <ul className="liste bloc-liste">
                  {l.map((e) => (
                    <li key={e.id}>
                      <LigneExercice exercice={e} etat={etatExercice(progression.exercices[e.id])} surExercice={surExercice} />
                    </li>
                  ))}
                </ul>
              </section>
            )
          })
        ) : (
          <p className="doux">Pas encore d’exercice pour cette notion.</p>
        ))}
    </main>
  )
}
