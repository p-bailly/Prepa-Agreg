// Mélange de Fisher-Yates avec un tirage injectable (Math.random par défaut, graine fixe dans les tests).
export function melanger<T>(t: readonly T[], alea: () => number = Math.random): T[] {
  const c = [...t]
  for (let i = c.length - 1; i > 0; i--) {
    const j = Math.floor(alea() * (i + 1))
    ;[c[i], c[j]] = [c[j], c[i]]
  }
  return c
}
