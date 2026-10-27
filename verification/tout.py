"""Lance tous les scripts de vérification calculatoire du contenu (SymPy, hors APK).

Chaque fichier verification/<theme>.py est fait d'assert commentés qui citent l'identifiant de la carte ou de
l'exercice vérifié. Voir CLAUDE.md, « Contrôle de l'exactitude ». Lancer : python3 verification/tout.py
"""
import pathlib
import runpy
import sys

dossier = pathlib.Path(__file__).parent
scripts = sorted(p for p in dossier.glob("*.py") if p.name != "tout.py")
for script in scripts:
    print(f"— {script.name}")
    runpy.run_path(str(script), run_name="__main__")
print(f"{len(scripts)} script(s) de vérification passés.")
sys.exit(0)
