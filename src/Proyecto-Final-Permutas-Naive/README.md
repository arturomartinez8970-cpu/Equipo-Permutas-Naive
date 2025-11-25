# Compilación del módulo `busqueda` (pybind11)

Instrucciones simples para compilar e instalar la extensión C++ (`busqueda`) usada por `main.py`.

Requisitos mínimos
- Windows con Python 3.x instalado.
- `pip` disponible.
- Visual Studio Build Tools (MSVC) con el workload "Desktop development with C++" instalado. Si no lo tiene, instale "Build Tools for Visual Studio".

Pasos (PowerShell)

1. Abra PowerShell en la carpeta del proyecto (la que contiene `setup.py`). Por ejemplo:

```powershell
cd "C:\Users\artur\Downloads\proyecto final_ Permutas Naive\proyecto final_ Permutas Naive"
```

2. Instale/actualice las dependencias de construcción:

```powershell
python -m pip install --upgrade pip setuptools wheel pybind11
```

3. Compile la extensión (comando fiable en PowerShell):

```powershell
python .\setup.py build_ext --inplace
```

Problemas comunes
- Error "'cl.exe' no se reconoce": no tiene MSVC en `PATH`. Abra el "Developer Command Prompt for VS" o instale Visual Studio Build Tools y luego compile desde ese prompt.
- Error relacionado con `pybind11` o includes: asegúrese de haber ejecutado `pip install pybind11` en el mismo intérprete de Python.
- Rutas con espacios: normalmente no es problema si está en la carpeta correcta; si falla, pruebe mover el proyecto a una ruta sin espacios.

Uso después de compilar
- Si la compilación fue correcta, en la carpeta aparecerá un archivo `busqueda*.pyd`. Puede ejecutar `main.py` que importa `busqueda` normalmente:

```powershell
python .\main.py
```

SSL y `certifi`
- Si al ejecutar `main.py` obtiene errores de certificado SSL, agregue este snippet al inicio de `main.py` (antes de usar `Entrez`):

```python
try:
		import certifi
		import os
		os.environ['SSL_CERT_FILE'] = certifi.where()
except Exception:
		# If certifi is not available yet, instale con:
		# python -m pip install certifi
		pass
```

- Para instalar `certifi`:

```powershell
python -m pip install certifi
```

Imagen (`adn.png`) — dónde cambiar la ruta
- El archivo `main.py` carga una imagen con `tk.PhotoImage(...)`. Busque la línea que asigna `img_adn`, por ejemplo:

```python
img_adn = tk.PhotoImage(file="C:\\Users\\migue\\Desktop\\proyectoFinal2\\adn.png").subsample(20,20)
```

- Opciones para hacer que el programa encuentre la imagen:
	- Opción 1 (recomendada): copie su `adn.png` dentro de la carpeta del proyecto (la misma que contiene `main.py`) y cambie la línea anterior a:

		```python
		img_adn = tk.PhotoImage(file="adn.png").subsample(20,20)
		```

	- Opción 2: sustituya la ruta absoluta por la ruta correcta donde está la imagen en su equipo. Use doble backslash `\\` en Windows o una cadena raw `r"C:\path\to\adn.png"`.

		```python
		img_adn = tk.PhotoImage(file=r"C:\ruta\correcta\a\adn.png").subsample(20,20)
		```

	- Opción 3: si no quiere usar imagen, puede dejar un placeholder (imagen vacía). El programa maneja un fallback si no existe la imagen.

Consejo rápido
- Si no está seguro de la ruta, abra el Explorador de Windows, vaya al archivo `adn.png`, haga click derecho → Propiedades → Copiar la ruta completa y péguela en el código usando `r"..."`.

