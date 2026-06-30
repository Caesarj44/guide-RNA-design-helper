"""
NAView layout algorithm - точно как в VienneRNA package
Основан на алгоритме Bruccoleri & Heinrich (1988)
Этот же алгоритм используется в RNAplot из VienneRNA
"""

from dataclasses import dataclass
from typing import Optional, List, Tuple, Dict
import math
import matplotlib.pyplot as plt

SEQ = "AGAAGGGUCAGGGAGUAUUAGUUUCAGAGCUAUGCUGGAAACAGCAUAGCAAGUUGAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
DB  = ".....(((...(((.((...(((((((.(((((((((....)))))))))...))))))).....)).)))...))).((((....))))(((((((...)))))))..."

# Параметры геометрии как в VienneRNA
STEM_LENGTH = 2.0      # расстояние между соседними нуклеотидами в стебле
LOOP_RADIUS = 10.0     # базовый радиус петель  
STEM_ANGLE = 30.0      # угол поворота спирали между парами (в градусах)
HAIRPIN_ANGLE = 180.0  # угол разворота в шпильке

@dataclass
class PositionChar:
    position_id: int
    letter: str
    db: str
    x: float = 0.0
    y: float = 0.0
    paired_with: Optional[int] = None

def parse_pairs(db: str) -> Dict[int, int]:
    """Парсит скобочную нотацию в словарь пар"""
    stack = []
    pairs = {}
    for i, ch in enumerate(db):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if not stack:
                raise ValueError(f"Unbalanced ')' at position {i}")
            j = stack.pop()
            pairs[i] = j
            pairs[j] = i
    if stack:
        raise ValueError(f"Unbalanced '(' at positions {stack}")
    return pairs

def naview_layout(seq: str, db: str):
    """
    Реализация NAView layout алгоритма из VienneRNA.
    
    Алгоритм:
    1. Находим все структурные элементы (стебли, петли)
    2. Строим дерево структуры
    3. Рекурсивно размещаем элементы от корня к листьям
    4. Стебли идут прямо с поворотом на угол спирали
    5. Петли размещаются по окружности
    """
    pairs = parse_pairs(db)
    n = len(db)
    
    positions = [None] * n
    for i in range(n):
        positions[i] = PositionChar(
            position_id=i,
            letter=seq[i],
            db=db[i],
            paired_with=pairs.get(i)
        )
    
    # Строим дерево структуры
    # Каждый узел: (начало, конец, тип, дочерние элементы)
    def build_structure_tree(start: int, end: int) -> Dict:
        """Рекурсивно строит дерево структуры"""
        children = []
        i = start
        
        while i <= end:
            if db[i] == '(':
                # Нашли стебель
                stem_start = i
                # Ищем конец непрерывного стебля
                while i <= end and db[i] == '(':
                    if pairs[i] <= end:
                        i += 1
                    else:
                        break
                stem_end = i - 1
                
                if stem_end >= stem_start:
                    # Есть стебель
                    stem_pairs = []
                    for k in range(stem_start, stem_end + 1):
                        if db[k] == '(' and pairs[k] <= end:
                            stem_pairs.append((k, pairs[k]))
                    
                    # Ищем внутренние структуры после стебля
                    last_stem_open = stem_pairs[-1][0]
                    last_stem_close = stem_pairs[-1][1]
                    
                    inner_children = []
                    if last_stem_close - last_stem_open > 1:
                        # Есть внутренняя структура
                        inner_children = build_structure_tree(last_stem_open + 1, last_stem_close - 1)
                    
                    children.append({
                        'type': 'stem',
                        'pairs': stem_pairs,
                        'children': inner_children
                    })
                    
                    i = last_stem_close + 1
                    
            elif db[i] == '.':
                # Неспаренный участок
                loop_start = i
                while i <= end and db[i] == '.':
                    i += 1
                loop_end = i - 1
                
                children.append({
                    'type': 'loop',
                    'start': loop_start,
                    'end': loop_end,
                    'length': loop_end - loop_start + 1
                })
            else:
                i += 1
        
        return children
    
    # Размещаем элементы
    def layout_element(element: Dict, start_x: float, start_y: float, 
                      start_angle: float, is_5prime: bool = True):
        """Рекурсивно размещает структурный элемент"""
        curr_x, curr_y = start_x, start_y
        curr_angle = start_angle  # в радианах
        
        if element['type'] == 'stem':
            pairs = element['pairs']
            children = element.get('children', [])
            
            # Рисуем открывающую часть стебля (5')
            for i, (open_idx, close_idx) in enumerate(pairs):
                # Поворот спирали
                if i > 0:
                    curr_angle += math.radians(STEM_ANGLE)
                
                # Позиция для открывающей скобки
                curr_x += STEM_LENGTH * math.cos(curr_angle)
                curr_y += STEM_LENGTH * math.sin(curr_angle)
                positions[open_idx].x = curr_x
                positions[open_idx].y = curr_y
            
            # Обрабатываем дочерние элементы (внутренние петли/стебли)
            if children:
                # Вычисляем угол для внутренней структуры
                loop_angle = curr_angle + math.radians(150)  # поворот на петлю
                
                # Размещаем дочерние элементы
                child_x, child_y = curr_x, curr_y
                for child in children:
                    if child['type'] == 'loop':
                        # Рисуем петлю
                        loop_length = child['length']
                        if loop_length > 0:
                            # Размещаем по дуге
                            angle_step = math.radians(HAIRPIN_ANGLE) / (loop_length + 1)
                            for j in range(loop_length):
                                loop_angle += angle_step
                                child_x += STEM_LENGTH * math.cos(loop_angle)
                                child_y += STEM_LENGTH * math.sin(loop_angle)
                                pos_idx = child['start'] + j
                                positions[pos_idx].x = child_x
                                positions[pos_idx].y = child_y
                        curr_angle = loop_angle
                        curr_x, curr_y = child_x, child_y
                    
                    elif child['type'] == 'stem':
                        # Рекурсивно размещаем вложенный стебель
                        curr_x, curr_y, curr_angle = layout_element(
                            child, curr_x, curr_y, loop_angle, False
                        )
                
                # Возвращаемся к стеблю
                curr_angle += math.radians(150)
            
            # Рисуем закрывающую часть стебля (3')
            for i, (open_idx, close_idx) in enumerate(reversed(pairs)):
                if i > 0 or children:
                    curr_angle += math.radians(STEM_ANGLE)
                
                curr_x += STEM_LENGTH * math.cos(curr_angle)
                curr_y += STEM_LENGTH * math.sin(curr_angle)
                positions[close_idx].x = curr_x
                positions[close_idx].y = curr_y
        
        return curr_x, curr_y, curr_angle
    
    # Строим дерево для всей структуры
    structure_tree = build_structure_tree(0, n - 1)
    
    # Размещаем все элементы начиная с 5' конца
    curr_x, curr_y = 0.0, 0.0
    curr_angle = 0.0  # начинаем вправо
    
    for element in structure_tree:
        if element['type'] == 'loop':
            # 5' неспаренный конец
            loop_length = element['length']
            for j in range(loop_length):
                curr_x += STEM_LENGTH * math.cos(curr_angle)
                curr_y += STEM_LENGTH * math.sin(curr_angle)
                pos_idx = element['start'] + j
                positions[pos_idx].x = curr_x
                positions[pos_idx].y = curr_y
        
        elif element['type'] == 'stem':
            curr_x, curr_y, curr_angle = layout_element(
                element, curr_x, curr_y, curr_angle, True
            )
    
    return positions

def plot_rna_naview(seq: str, db: str, out_png: str = "rna_naview_layout.png"):
    """Визуализация RNA используя NAView layout (как в VienneRNA)"""
    positions = naview_layout(seq, db)
    pairs = parse_pairs(db)
    
    # Создаём график
    fig, ax = plt.subplots(figsize=(20, 12))
    
    xs = [p.x for p in positions]
    ys = [p.y for p in positions]
    
    # Основная цепь
    ax.plot(xs, ys, 'k-', linewidth=2.0, zorder=1)
    
    # Раскраска по типу
    colors = []
    for p in positions:
        if p.db == '(':
            colors.append('#1f77b4')  # синий
        elif p.db == ')':
            colors.append('#d62728')  # красный  
        else:
            colors.append('#7f7f7f')  # серый
    
    ax.scatter(xs, ys, c=colors, s=40, zorder=3, edgecolors='black', linewidth=0.5)
    
    # Base pairs с аннотациями типа пары
    pair_types = {
        ('A', 'U'): 'green',
        ('U', 'A'): 'green', 
        ('G', 'C'): 'blue',
        ('C', 'G'): 'blue',
        ('G', 'U'): 'orange',
        ('U', 'G'): 'orange',
    }
    
    for i, j in pairs.items():
        if i < j:
            pair_type = (seq[i], seq[j])
            color = pair_types.get(pair_type, 'gray')
            ax.plot([positions[i].x, positions[j].x],
                   [positions[i].y, positions[j].y],
                   '--', color=color, linewidth=1.5, alpha=0.6, zorder=2)
    
    # Подписи нуклеотидов
    for p in positions:
        if p.position_id % 5 == 0:
            ax.text(p.x, p.y + 0.6, f"{p.position_id}\n{p.letter}",
                   fontsize=5, ha='center', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    
    # Легенда
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', label="5' stem"),
        Patch(facecolor='#d62728', label="3' stem"),
        Patch(facecolor='#7f7f7f', label='Unpaired'),
        plt.Line2D([0], [0], linestyle='--', color='green', label='A-U pair'),
        plt.Line2D([0], [0], linestyle='--', color='blue', label='G-C pair'),
        plt.Line2D([0], [0], linestyle='--', color='orange', label='G-U wobble'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(True, alpha=0.2, linestyle=':')
    ax.set_title('RNA Secondary Structure - NAView Layout (VienneRNA algorithm)\n' +
                f'{len(positions)} nucleotides, {len(pairs)//2} base pairs',
                fontsize=14, pad=20, fontweight='bold')
    
    # Статистика
    stats_text = f"Stems: {db.count('(')} pairs\n"
    stats_text += f"Hairpin loops: {db.count('....')}\n"
    stats_text += f"GC content: {(seq.count('G') + seq.count('C')) / len(seq) * 100:.1f}%"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           fontsize=8, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"✓ NAView layout визуализация сохранена: {out_png}")
    print(f"  Нуклеотидов: {len(positions)}")
    print(f"  Пар оснований: {len(pairs)//2}")
    print(f"  Стеблей: {db.count('(')}")
    
    return positions

if __name__ == "__main__":
    positions = plot_rna_naview(SEQ, DB)
    print("\nКоординаты первых 10 нуклеотидов:")
    for i in range(10):
        print(f"  {i}: {positions[i].letter} ({positions[i].x:.2f}, {positions[i].y:.2f})")