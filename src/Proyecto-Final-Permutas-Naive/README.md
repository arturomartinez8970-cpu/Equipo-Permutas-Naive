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

Uso después de compilar
- Si la compilación fue correcta, en la carpeta aparecerá un archivo `busqueda*.pyd`. Puede ejecutar `main.py` que importa `busqueda` normalmente:

```powershell
python .\main.py
```

### Error de certificacion
SSL y `certifi`
- Si al ejecutar `main.py` obtiene errores de certificado SSL, agregue este snippet al inicio de `main.py` (luego de `import os`):

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

### Problema de imagen (`adn.png`) 
— dónde cambiar la ruta
- El archivo `main.py` carga una imagen con `tk.PhotoImage(...)`. Busque la línea que asigna `img_adn`, por ejemplo:

```python
img_adn = tk.PhotoImage(file="C:\\Users\\migue\\Desktop\\proyectoFinal2\\adn.png").subsample(20,20)
```
- Como hacer que el programa encuentre la imagen:
	- Como adn ya esta en la carpeta del proyecto (la misma que contiene `main.py`) solo habra que cambiar la línea anterior a:

		```python
		img_adn = tk.PhotoImage(file="adn.png").subsample(20,20)
		```

