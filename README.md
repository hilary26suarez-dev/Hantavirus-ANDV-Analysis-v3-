# Hantavirus ANDV: Computational Mutation Analysis 🦠

### 🧬 Resumen del Proyecto (Overview)
Esta es una iniciativa de investigación personal nacida de la curiosidad científica, centrada en el análisis computacional de mutaciones puntuales en las glicoproteínas (Gn/Gc) del **Virus Andes (ANDV)**. El objetivo es evaluar cómo sustituciones específicas de aminoácidos alteran las propiedades electrostáticas del virus y su afinidad potencial por el receptor humano **Integrina β₃**.

### 🧪 Hipótesis
Sustituciones que incrementan la carga neta positiva en la proteína Gn (ej. **E50R**) mejoran la afinidad de unión a la integrina β₃ humana (que posee carga negativa), facilitando potencialmente la entrada viral y aumentando la eficiencia de infectividad.

### 💻 Características Técnicas (Features)
* **Recuperación Automatizada de Datos:** Uso de `Biopython` y la API `Entrez` para obtener secuencias de proteínas directamente del NCBI.
* **Mutagénesis In Silico:** Simulación de sustituciones específicas en posiciones clave (50, 100, 150).
* **Perfilado Fisicoquímico:** Cálculo de carga neta, densidad de carga y propensión a estructuras secundarias (Hélice/Lámina-β).
* **Análisis de Hidropatía:** Evaluación del potencial de contacto superficial basado en la escala de Kyte-Doolittle.

### 📊 Metodología
1.  **Adquisición de Secuencias:** Análisis de UniProt P05106 (Integrina β₃ humana) y glicoproteínas de ANDV.
2.  **Simulación Computacional:** Análisis basado en Python de 9 mutaciones candidatas.
3.  **Visualización de Datos:** Generación de datasets en JSON y reportes interactivos en HTML.

### 🛡️ Declaración Ética y de Bioseguridad
Esta investigación se realiza con fines **exclusivamente preventivos, defensivos y educativos**. El objetivo es profundizar en el conocimiento de la patogenicidad viral para apoyar el futuro desarrollo de terapias antivirales y herramientas de diagnóstico. 
* **No implica manipulación física de patógenos.**
* Se basa en datos de secuencias de acceso público.
* Busca identificar blancos terapéuticos para el control de brotes.

---

### 🛠️ Estructura del Repositorio
* `hantavirus_analysis.py`: Script principal de análisis bioinformático.
* `datos_20260508_234044.json`: Dataset con los resultados de las mutaciones simuladas.
* `reporte_20260508_234044.html`: Reporte visual de hallazgos.
* `REPORTE_FINAL_ANÁLISIS_HANTAVIRUS.txt`: Documentación detallada del proyecto.

---
**Autor:** Hilary Gretchan Suárez Fonseca
**Estado:** Iniciativa de Investigación Independiente
**Institución de Referencia:** UCIMED
