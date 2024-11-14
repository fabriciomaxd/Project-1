#!/usr/bin/env python3
# -- coding: UTF-8 --

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: <SEU NOME AQUI>
Disciplina: Computação Gráfica
Data: <DATA DE INÍCIO DA IMPLEMENTAÇÃO>
"""

import time         # Para operações com tempo
import gpu          # Simula os recursos de uma GPU
import math         # Funções matemáticas
import numpy as np  # Biblioteca do Numpy

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela  
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante
    projection_matrix = np.identity(4)
    model_matrix = np.identity(4)
    view_matrix = np.identity(4)

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far

    @staticmethod
    def polypoint2D(point, colors):
        """Função usada para renderizar Polypoint2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polypoint2D
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
        # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
        # você pode assumir inicialmente o desenho dos pontos com a cor emissiva (emissiveColor).


        emissive_color = colors.get('emissiveColor')  # Usa branco como padrão
        emissive_color = [int(c * 255) for c in emissive_color]

        for i in range(0, len(point), 2):
            x = int(point[i])        # Coordenada x
            y = int(point[i + 1])    # Coordenada y
            gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8,emissive_color)


    @staticmethod
    def polyline2D(lineSegments, colors):
        """Função usada para renderizar Polyline2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polyline2D
        # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
        # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
        # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
        # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
        # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
        # vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        emissive_color = colors.get('emissiveColor') 
        emissive_color = [int(c * 255) for c in emissive_color]

        def bres(x0, y0, x1, y1):
            pontos = []
            dx = abs(x1 - x0)
            dy = abs(y1 - y0)
            sx = 1 if x0 < x1 else -1
            sy = 1 if y0 < y1 else -1
            erro = dx - dy

            while True:
                pontos.append((x0, y0))  
                if x0 == x1 and y0 == y1:
                    break
                erro2 = erro * 2
                if erro2 > -dy:
                    erro -= dy
                    x0 += sx
                if erro2 < dx:
                    erro += dx
                    y0 += sy

            return pontos

        for i in range(0, len(lineSegments) - 2, 2):
                x0 = int(lineSegments[i])        
                y0 = int(lineSegments[i + 1])    
                x1 = int(lineSegments[i + 2])   
                y1 = int(lineSegments[i + 3])    
                
                line_points = bres(x0, y0, x1, y1)
                
                for point in line_points:
                    gpu.GPU.draw_pixel(point, gpu.GPU.RGB8, emissive_color)

    @staticmethod
    def circle2D(radius, colors):
        """Função usada para renderizar Circle2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Circle2D
        # Nessa função você receberá um valor de raio e deverá desenhar o contorno de
        # um círculo.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Circle2D
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

            
        emissive_color = colors.get('emissiveColor', [1, 0, 0])  # Default para vermelho
        emissive_color = [int(c * 255) for c in emissive_color]  # Converte para [0, 255]

        # Meio da Tela
        pos_x = GL.width // 2
        pos_y = GL.height // 2

        # Verifique se o raio é válido
        if radius <= 0:
            raise ValueError("Radius must be greater than zero.")

        # Variáveis
        x = radius
        y = 0
        p = 1 - radius  

        # Dimensões do framebuffer
        fb_width, fb_height = GL.width, GL.height

        # Desenhar o círculo
        while x >= y:
            # Verificar e desenhar os pontos simétricos
            if 0 <= pos_x + x < fb_width and 0 <= pos_y + y < fb_height:
                gpu.GPU.draw_pixel([int(pos_x + x), int(pos_y + y)], gpu.GPU.RGB8, emissive_color)  # Quadrante I
            if 0 <= pos_x - x < fb_width and 0 <= pos_y + y < fb_height:
                gpu.GPU.draw_pixel([int(pos_x - x), int(pos_y + y)], gpu.GPU.RGB8, emissive_color)  # Quadrante II
            if 0 <= pos_x + x < fb_width and 0 <= pos_y - y < fb_height:
                gpu.GPU.draw_pixel([int(pos_x + x), int(pos_y - y)], gpu.GPU.RGB8, emissive_color)  # Quadrante III
            if 0 <= pos_x - x < fb_width and 0 <= pos_y - y < fb_height:
                gpu.GPU.draw_pixel([int(pos_x - x), int(pos_y - y)], gpu.GPU.RGB8, emissive_color)  # Quadrante IV
            if 0 <= pos_x + y < fb_width and 0 <= pos_y + x < fb_height:
                gpu.GPU.draw_pixel([int(pos_x + y), int(pos_y + x)], gpu.GPU.RGB8, emissive_color)  # Quadrante I
            if 0 <= pos_x - y < fb_width and 0 <= pos_y + x < fb_height:
                gpu.GPU.draw_pixel([int(pos_x - y), int(pos_y + x)], gpu.GPU.RGB8, emissive_color)  # Quadrante II
            if 0 <= pos_x + y < fb_width and 0 <= pos_y - x < fb_height:
                gpu.GPU.draw_pixel([int(pos_x + y), int(pos_y - x)], gpu.GPU.RGB8, emissive_color)  # Quadrante III
            if 0 <= pos_x - y < fb_width and 0 <= pos_y - x < fb_height:
                gpu.GPU.draw_pixel([int(pos_x - y), int(pos_y - x)], gpu.GPU.RGB8, emissive_color)  # Quadrante IV

            y += 1
            
            # Atualizar o valor de p
            if p <= 0:
                p += 2 * y + 1  # Move para o leste
            else:
                x -= 1  # Move para o sudoeste
                p += 2 * y - 2 * x + 1  # Move para o leste


    @staticmethod
    def triangleSet2D(vertices, colors):
            """Função usada para renderizar TriangleSet2D."""
            # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#TriangleSet2D
            # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
            # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
            # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
            # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
            # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
            # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
            # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).
            # Itera sobre os pontos em conjuntos de 3, pois cada triângulo precisa de 3 pontos (6 valores no total)

            # Cor emissiva para o triângulo, convertida para formato 0-255
            emissive_color = colors.get('emissiveColor', (1.0, 1.0, 1.0))  # Padrão é branco
            emissive_color = [int(c * 255) for c in emissive_color]

            # Itera sobre os vértices em grupos de 6 para cada triângulo
            for i in range(0, len(vertices), 6):
                # Extrai as coordenadas dos vértices do triângulo, convertendo para inteiros
                x0, y0 = int(vertices[i]), int(vertices[i + 1])
                x1, y1 = int(vertices[i + 2]), int(vertices[i + 3])
                x2, y2 = int(vertices[i + 4]), int(vertices[i + 5])

                # Calcula a área total do triângulo usando a fórmula da área
                area = 0.5 * abs(x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1))

                # Define a bounding box para o triângulo
                min_x = max(min(x0, x1, x2), 0)
                max_x = min(max(x0, x1, x2), GL.width - 1)
                min_y = max(min(y0, y1, y2), 0)
                max_y = min(max(y0, y1, y2), GL.height - 1)

                # Itera sobre os pixels dentro da bounding box
                for x in range(min_x, max_x + 1):
                    for y in range(min_y, max_y + 1):
                        # Calcula as áreas dos sub-triângulos
                        area1 = 0.5 * abs(x0 * (y1 - y) + x1 * (y - y0) + x * (y0 - y1))
                        area2 = 0.5 * abs(x1 * (y2 - y) + x2 * (y - y1) + x * (y1 - y2))
                        area3 = 0.5 * abs(x2 * (y0 - y) + x0 * (y - y2) + x * (y2 - y0))

                        # Verifica se a soma das áreas dos sub-triângulos é aproximadamente igual à área do triângulo
                        if abs(area - (area1 + area2 + area3)) < 1e-6:
                            # Se o pixel está dentro do triângulo, define sua cor
                            gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, emissive_color)
                
    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleSet
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, você pode assumir
        # inicialmente, para o TriangleSet, o desenho das linhas com a cor emissiva
        # (emissiveColor), conforme implementar novos materias você deverá suportar outros
        # tipos de cores.

        
        # Cor
        #emissive_color = colors.get('emissiveColor', (1.0, 1.0, 1.0))
        #emissive_color = [int(c * 255) for c in emissive_color[:3]]pyt
        #emissive_color = [max(0, min(255, c)) for c in emissive_color]

        #Pontos do triangulo 
        for i in range(0, len(point), 9):
            x0, y0, z0 = point[i], point[i + 1], point[i + 2]
            x1, y1, z1 = point[i + 3], point[i + 4], point[i + 5]
            x2, y2, z2 = point[i + 6], point[i + 7], point[i + 8]

            # Projeção
            p0 = np.dot(GL.projection_matrix, np.dot(GL.view_matrix, np.dot(GL.model_matrix, np.array([x0, y0, z0, 1]))))
            p1 = np.dot(GL.projection_matrix, np.dot(GL.view_matrix, np.dot(GL.model_matrix, np.array([x1, y1, z1, 1]))))
            p2 = np.dot(GL.projection_matrix, np.dot(GL.view_matrix, np.dot(GL.model_matrix, np.array([x2, y2, z2, 1]))))

            # Coordenadas 2D
            p0 /= p0[3]
            p1 /= p1[3]
            p2 /= p2[3]

            # Coordenadas tela
            vertices_2d = [
                int((p0[0] + 1) * GL.width / 2), int((1 - p0[1]) * GL.height / 2),
                int((p1[0] + 1) * GL.width / 2), int((1 - p1[1]) * GL.height / 2),
                int((p2[0] + 1) * GL.width / 2), int((1 - p2[1]) * GL.height / 2)
            ]

            # Limites
            vertices_2d = [
                max(0, min(GL.width - 1, vertices_2d[0])), max(0, min(GL.height - 1, vertices_2d[1])),
                max(0, min(GL.width - 1, vertices_2d[2])), max(0, min(GL.height - 1, vertices_2d[3])),
                max(0, min(GL.width - 1, vertices_2d[4])), max(0, min(GL.height - 1, vertices_2d[5]))
            ]

            
            #print(emissive_color)
            # Triangulos
            #vertices_2d = [1, 0, 0, 1, 0, 1]
            #print(vertices_2d)
            GL.triangleSet2D(vertices_2d,colors)
            

    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
              # Matriz de projeção perspectiva
        # Cria a matriz de projeção perspectiva
      
        # Cria a matriz de projeção perspectiva
        # Configura a matriz de projeção perspectiva

    # Matriz de projeção perspectiva
        aspect_ratio = GL.width / GL.height
        near = 0.1
        far = 1000.0
        f = 1 / np.tan(fieldOfView / 2) 
        
        GL.projection_matrix = np.array([
            [f / aspect_ratio, 0, 0, 0],
            [0, f, 0, 0],
            [0, 0, (far + near) / (near - far), (2 * far * near) / (near - far)],
            [0, 0, -1, 0]
        ])

        # Array Numpy
        position = np.array(position)
        tx, ty, tz = -position

        # Componentes 
        ux, uy, uz, angle = orientation
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        one_minus_cos = 1 - cos_a

        # Rotação
        rotation_matrix = np.array([
            [cos_a + ux**2 * one_minus_cos, ux * uy * one_minus_cos - uz * sin_a, ux * uz * one_minus_cos + uy * sin_a, 0],
            [uy * ux * one_minus_cos + uz * sin_a, cos_a + uy**2 * one_minus_cos, uy * uz * one_minus_cos - ux * sin_a, 0],
            [uz * ux * one_minus_cos - uy * sin_a, uz * uy * one_minus_cos + ux * sin_a, cos_a + uz**2 * one_minus_cos, 0],
            [0, 0, 0, 1]
        ])

        # Translação
        translation_matrix = np.array([
            [1, 0, 0, tx],
            [0, 1, 0, ty],
            [0, 0, 1, tz],
            [0, 0, 0, 1]
        ])

        # Matriz de visualização final
        GL.view_matrix = np.dot(rotation_matrix, translation_matrix)

    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo para depois potencialmente usar em outras chamadas. 
        # Quando começar a usar Transforms dentre de outros Transforms, mais a frente no curso
        # Você precisará usar alguma estrutura de dados pilha para organizar as matrizes.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.

        # Matriz de escala
        scale_matrix = np.array([
            [scale[0], 0, 0, 0],
            [0, scale[1], 0, 0],
            [0, 0, scale[2], 0],
            [0, 0, 0, 1]
        ])
        
        # Matriz de translação
        translation_matrix = np.array([
            [1, 0, 0, translation[0]],
            [0, 1, 0, translation[1]],
            [0, 0, 1, translation[2]],
            [0, 0, 0, 1]
        ])

        # Matriz de rotação com base no eixo e ângulo fornecidos
        ux, uy, uz, angle = rotation
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        one_minus_cos = 1 - cos_a

        rotation_matrix = np.array([
            [cos_a + ux**2 * one_minus_cos, ux * uy * one_minus_cos - uz * sin_a, ux * uz * one_minus_cos + uy * sin_a, 0],
            [uy * ux * one_minus_cos + uz * sin_a, cos_a + uy**2 * one_minus_cos, uy * uz * one_minus_cos - ux * sin_a, 0],
            [uz * ux * one_minus_cos - uy * sin_a, uz * uy * one_minus_cos + ux * sin_a, cos_a + uz**2 * one_minus_cos, 0],
            [0, 0, 0, 1]
        ])

        # Composição da matriz de transformação final do objeto
        GL.model_matrix = np.dot(translation_matrix, np.dot(rotation_matrix, scale_matrix))

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Saindo de Transform")

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleStripSet
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        #emissive_color = colors.get('emissiveColor', (1.0, 1.0, 1.0))
        #emissive_color = [int(c * 255) for c in emissive_color]

        #emissive_color = colors.get('emissiveColor', (1.0, 1.0, 1.0))
        #emissive_color = [int(c * 255) for c in emissive_color]

        offset = 0
        for count in stripCount:
            for i in range(count - 2):
                idx1, idx2, idx3 = offset + i, offset + i + 1, offset + i + 2
                vertices = [
                    point[idx1 * 3], point[idx1 * 3 + 1], point[idx1 * 3 + 2],
                    point[idx2 * 3], point[idx2 * 3 + 1], point[idx2 * 3 + 2],
                    point[idx3 * 3], point[idx3 * 3 + 1], point[idx3 * 3 + 2]
                ]
                # Use triangleSet para renderizar em 3D
                GL.triangleSet(vertices, colors)
            offset += count

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#IndexedTriangleStripSet
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        #emissive_color = colors.get('emissiveColor', (1.0, 1.0, 1.0))
        #emissive_color = [int(c * 255) for c in emissive_color]

        strip = []
        for i in index:
            if i == -1:
                # Conecta os pontos da faixa corrente
                for j in range(len(strip) - 2):
                    idx1, idx2, idx3 = strip[j], strip[j + 1], strip[j + 2]
                    vertices = [
                        point[idx1 * 3], point[idx1 * 3 + 1], point[idx1 * 3 + 2],
                        point[idx2 * 3], point[idx2 * 3 + 1], point[idx2 * 3 + 2],
                        point[idx3 * 3], point[idx3 * 3 + 1], point[idx3 * 3 + 2]
                    ]
                    # Use triangleSet para renderizar em 3D
                    if len(vertices) == 9:
                        GL.triangleSet(vertices, colors)
                strip = []  # Reseta para a próxima faixa
            else:
                strip.append(i)
    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#IndexedFaceSet
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão não possui uma ordem oficial, mas em geral se o primeiro ponto com os dois
        # seguintes e depois este mesmo primeiro ponto com o terçeiro e quarto ponto. Por exemplo: numa
        # sequencia 0, 1, 2, 3, 4, -1 o primeiro triângulo será com os vértices 0, 1 e 2, depois serão
        # os vértices 0, 2 e 3, e depois 0, 3 e 4, e assim por diante, até chegar no final da lista.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        """Função usada para renderizar IndexedFaceSet."""
        #emissive_color = colors.get('emissiveColor', (1.0, 1.0, 1.0))
        #emissive_color = [int(c * 255) for c in emissive_color]

        face = []
        for i in coordIndex:
            if i == -1:
                # Conecta os pontos da face atual em triângulos
                for j in range(1, len(face) - 1):
                    idx1, idx2, idx3 = face[0], face[j], face[j + 1]
                    vertices = [
                        coord[idx1 * 3], coord[idx1 * 3 + 1], coord[idx1 * 3 + 2],
                        coord[idx2 * 3], coord[idx2 * 3 + 1], coord[idx2 * 3 + 2],
                        coord[idx3 * 3], coord[idx3 * 3 + 1], coord[idx3 * 3 + 2]
                    ]
                    # Use triangleSet para renderizar em 3D
                    GL.triangleSet(vertices, colors)
                face = []  # Reseta para a próxima face
            else:
                face.append(i)

    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone
        # A função cone é usada para desenhar cones na cena. O cone é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento bottomRadius especifica o
        # raio da base do cone e o argumento height especifica a altura do cone.
        # O cone é alinhado com o eixo Y local. O cone é fechado por padrão na base.
        # Para desenha esse cone você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cone : bottomRadius = {0}".format(bottomRadius)) # imprime no terminal o raio da base do cone
        print("Cone : height = {0}".format(height)) # imprime no terminal a altura do cone
        print("Cone : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder
        # A função cylinder é usada para desenhar cilindros na cena. O cilindro é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da base do cilindro e o argumento height especifica a altura do cilindro.
        # O cilindro é alinhado com o eixo Y local. O cilindro é fechado por padrão em ambas as extremidades.
        # Para desenha esse cilindro você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cylinder : radius = {0}".format(radius)) # imprime no terminal o raio do cilindro
        print("Cylinder : height = {0}".format(height)) # imprime no terminal a altura do cilindro
        print("Cylinder : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#DirectionalLight
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
