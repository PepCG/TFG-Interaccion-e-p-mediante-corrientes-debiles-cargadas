# TFG: Interaccion e-p mediante corrientes debiles cargadas
Se incluyen todos los códigos en Python asociados al TFG, a partir de los cuales es posible reproducir los resultados obtenidos. 

## Código principal
En primer lugar, se presenta el código principal, en un archivo llamado "cross_section_nucleon.py". Dicho código nos permite calcular y plotear las secciones eficaces (en función del ángulo de dispersión, energía transferida y cuadrimomento al cuadrado transferido en la interacción) de los procesos descritos en el TFG. Para su uso, este necesita el archivo llamado "Input.txt", que nos permite definir los diferentes parámetros del proceso de dispersión que queramos caracterizar (para mayor info, leer "README_INPUT.md"). 

## Carpeta Resultados
En la carpeta llamada "resultados", se exponen los diferentes códigos que nos permiten obtener las gráficas mostradas en la Sección 6 del TFG. Cada código está contenido en una carpeta cuyo nombre hace referencia a la Figura que nos permitiría dibujar (es decir, la carpeta llamada "figura_6.1" nos permite graficar dicha figura). En este caso, será necesario guardar los datos a partir del código principal "cross_section_nucleon.py", haciendo uso de las líneas 540-546 comentadas en el mismo (se ha decicido no compartir los datos para dar una mayor maniobrabilidad al usuario y para que el repositorio quede más limpio y comprensible). A partir de dichos datos, siempre que sean almacenados en la carpeta donde se encuentra el código y nombrados en la forma inidicada, será posible graficar todas las figuras de la Sección 6 del TFG.
