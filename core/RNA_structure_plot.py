from __future__ import annotations

import math
from dataclasses import dataclass, field

import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans','DejaVu Sans']
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


PAIR_LEN = 1.3         
BASE_LEN = 1.3         
HELIX_STEP = 1.3    
MIN_LOOP_RADIUS = 1.8 
UNPAIRED_WEIGHT = 1.0   
SUBTREE_WEIGHT_EXP = 0.3  

BASE_COLORS = {
    'A': '#6fbf73',
    'U': '#f2c46d',
    'T': '#f2c46d',
    'G': '#e57373',
    'C': '#64b5f6',
}

def parse_dot_bracket(structure: str) -> list[int]:
    n = len(structure)
    pairs = [-1] * n
    bracket_pairs = {')': '(', ']': '[', '}': '{', '>': '<'}
    stacks: dict[str, list[int]] = {b: [] for b in '([{<'}

    for i, ch in enumerate(structure):
        if ch in stacks:
            stacks[ch].append(i)
        elif ch in bracket_pairs:
            opening = bracket_pairs[ch]
            if not stacks[opening]:
                raise ValueError(f"Unbalanced '{ch}' at position {i}")
            j = stacks[opening].pop()
            pairs[i] = j
            pairs[j] = i
        elif ch == '.':
            pass
        else:
            raise ValueError(f"Unknown symbol '{ch}' at position {i}")

    for b, stack in stacks.items():
        if stack:
            raise ValueError(f"Unclosed '{b}' at position {stack[0]}")

    return pairs



@dataclass
class Helix:
    """Stem: list of (i, j) pairs, i increasing, j decreasing."""
    base_pairs: list[tuple[int, int]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.base_pairs)


@dataclass
class Loop:
    unpaired: list[int] = field(default_factory=list)
    children: list[tuple[Helix, "Loop"]] = field(default_factory=list)
    order: list[tuple] = field(default_factory=list)
    is_exterior: bool = False


def build_tree(pairs: list[int]) -> Loop:
    n = len(pairs)
    root = Loop(is_exterior=True)

    def parse_loop(start: int, end: int, loop: Loop) -> None:
        i = start
        while i < end:
            j = pairs[i]
            if j == -1:
                loop.unpaired.append(i)
                loop.order.append(('nt', i))
                i += 1
            else:
                helix = Helix()
                a, b = i, j
                while True:
                    helix.base_pairs.append((a, b))
                    na, nb = a + 1, b - 1
                    if na < nb and pairs[na] == nb:
                        a, b = na, nb
                        continue
                    break
                child_loop = Loop()
                loop.children.append((helix, child_loop))
                loop.order.append(('helix', len(loop.children) - 1))
                parse_loop(a + 1, b, child_loop)
                i = j + 1

    parse_loop(0, n, root)
    return root


@dataclass
class Layout:
    coords: dict[int, tuple[float, float]] = field(default_factory=dict)
    pairs: list[int] = field(default_factory=list)


def loop_perimeter_units(loop: Loop) -> int:
    return len(loop.unpaired) + len(loop.children)


def _child_subtree_weight(helix: Helix, child: Loop) -> int:
    w = 2 * helix.length + len(child.unpaired)
    for h2, c2 in child.children:
        w += _child_subtree_weight(h2, c2)
    return w


def loop_radius(loop: Loop) -> float:
    units = loop_perimeter_units(loop)
    if units <= 2:
        return MIN_LOOP_RADIUS

    weighted_units = len(loop.unpaired) * UNPAIRED_WEIGHT
    for helix, child in loop.children:
        w = _child_subtree_weight(helix, child)
        weighted_units += max(1.5, w ** SUBTREE_WEIGHT_EXP)

    if weighted_units <= 2:
        return MIN_LOOP_RADIUS

    r = (BASE_LEN / 2) / math.sin(math.pi / weighted_units)
    return max(r * 1.15, MIN_LOOP_RADIUS)


def compute_layout(structure: str) -> Layout:
    pairs = parse_dot_bracket(structure)
    root = build_tree(pairs)
    coords: dict[int, tuple[float, float]] = {}

    def subtree_weight(loop: Loop) -> int:
        w = len(loop.unpaired)
        for helix, child in loop.children:
            w += 2 * helix.length + subtree_weight(child)
        return w

    def place_helix_chain(
            helix: Helix, start: tuple[float, float],
            direction: tuple[float, float],
            i_side: float
            ) -> tuple[float, float]:

        dx, dy = direction
        perp_x, perp_y = -dy, dx
        half = PAIR_LEN / 2
        for k, (i, j) in enumerate(helix.base_pairs):
            cx = start[0] + dx * (k * HELIX_STEP)
            cy = start[1] + dy * (k * HELIX_STEP)
            coords[i] = (cx + perp_x * half * i_side, cy + perp_y * half * i_side)
            coords[j] = (cx - perp_x * half * i_side, cy - perp_y * half * i_side)
        n = helix.length
        return (start[0] + dx * (n * HELIX_STEP),
                start[1] + dy * (n * HELIX_STEP))

    def place_internal_loop(loop: Loop, center: tuple[float, float],
                            entry_angle: float) -> None:

        r = loop_radius(loop)
        half_pair_angle = math.asin(min(1.0, (PAIR_LEN / 2) / r))
        remaining_arc = 2 * math.pi - 2 * half_pair_angle

        token_weights = []
        for token in loop.order:
            if token[0] == 'nt':
                token_weights.append(UNPAIRED_WEIGHT)
            else:
                helix, child = loop.children[token[1]]
                w = _child_subtree_weight(helix, child)
                token_weights.append(max(1.5, w ** SUBTREE_WEIGHT_EXP))

        total_w = sum(token_weights) if token_weights else 1.0
        margin_w = total_w / max(len(token_weights), 1) if token_weights else total_w
        unit_angle = remaining_arc / (total_w + margin_w) if (total_w + margin_w) > 0 else 0.0

        angle = entry_angle + half_pair_angle + (margin_w / 2) * unit_angle
        for idx, token in enumerate(loop.order):
            angle_center = angle + (token_weights[idx] / 2) * unit_angle
            x = center[0] + r * math.cos(angle_center)
            y = center[1] + r * math.sin(angle_center)
            if token[0] == 'nt':
                coords[token[1]] = (x, y)
            else:
                helix, child_loop = loop.children[token[1]]
                direction = (math.cos(angle_center), math.sin(angle_center))
                base_point = (x, y)
                end_point = place_helix_chain(helix, base_point, direction, i_side=-1.0)
                child_r = loop_radius(child_loop)*0.25
                child_center = (end_point[0] + direction[0] * child_r,
                                end_point[1] + direction[1] * child_r)
                place_internal_loop(child_loop, child_center, angle_center + math.pi)
            angle += token_weights[idx] * unit_angle

    def place_exterior_circular(loop: Loop, center: tuple[float, float] = (0.0, 0.0)) -> None:
        if len(loop.order) == 0:
            return

        r = loop_radius(loop)

        token_weights = []
        for token in loop.order:
            if token[0] == 'nt':
                token_weights.append(UNPAIRED_WEIGHT)
            else:
                helix, child = loop.children[token[1]]
                w = _child_subtree_weight(helix, child)
                token_weights.append(max(1.5, w ** SUBTREE_WEIGHT_EXP))

        total_w = sum(token_weights) if token_weights else 1.0
        margin_w = total_w / max(len(token_weights), 1) if token_weights else 1.0
        unit_angle = 2 * math.pi / (total_w + margin_w)

        angle = -math.pi / 2 + (margin_w / 2) * unit_angle
        for idx, token in enumerate(loop.order):
            angle_center = angle + (token_weights[idx] / 2) * unit_angle
            x = center[0] + r * math.cos(angle_center)
            y = center[1] + r * math.sin(angle_center)
            if token[0] == 'nt':
                coords[token[1]] = (x, y)
            else:
                helix, child_loop = loop.children[token[1]]
                direction = (math.cos(angle_center), math.sin(angle_center))
                base_point = (x, y)
                end_point = place_helix_chain(helix, base_point, direction, i_side=-1.0)
                child_r = loop_radius(child_loop)*0.25
                child_center = (end_point[0] + direction[0] * child_r,
                                end_point[1] + direction[1] * child_r)
                place_internal_loop(child_loop, child_center, angle_center + math.pi)
            angle += token_weights[idx] * unit_angle

    place_exterior_circular(root)

    for i in coords:
        x, y = coords[i]
        coords[i] = (-x, -y)

    return Layout(coords=coords, pairs=pairs)


def plot_rna_structure(
    sequence: str,
    structure: str,
    ax: plt.Axes | None = None,
    *,
    node_radius: float = 0.5,
    font_size: int = 8,
    show_indices: bool = False,
    title: str | None = None,
    figsize: tuple = (10, 10),
):
    if len(sequence) != len(structure):
        raise ValueError(
            f"Sequence length ({len(sequence)}) != structure length ({len(structure)})"
        )

    layout = compute_layout(structure)
    coords = layout.coords
    pairs = layout.pairs
    n = len(sequence)

    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor='none')
        created_fig = True
    else:
        fig = ax.figure

    xs = [coords[i][0] for i in range(n)]
    ys = [coords[i][1] for i in range(n)]

    ax.plot(xs, ys, '-', color='#888888', linewidth=1.0, zorder=1)

    drawn = set()
    for i, j in enumerate(pairs):
        if j != -1 and i < j and (i, j) not in drawn:
            drawn.add((i, j))
            x1, y1 = coords[i]
            x2, y2 = coords[j]
            ax.plot([x1, x2], [y1, y2], '-', color='#bbbbbb',
                    linewidth=1.0, zorder=1)

    for i in range(n):
        x, y = coords[i]
        base = sequence[i].upper()
        color = BASE_COLORS.get(base, '#cccccc')
        circle = Circle((x, y), node_radius, facecolor=color,
                        edgecolor='#333333', linewidth=1.0, zorder=2)
        ax.add_patch(circle)
        ax.text(x, y, base, ha='center', va='center',
                fontsize=font_size, zorder=3, fontweight='bold',
                color='#222222')
        if show_indices and (i % 10 == 0 or i == 0 or i == n - 1):
            ax.text(x, y - node_radius - 0.45, str(i + 1),
                    ha='center', va='center', fontsize=font_size * 0.7,
                    color='#555555', zorder=3)

    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title)

    if created_fig:
        fig.tight_layout()

    return fig, ax


def plot_rna(
    sequence: str,
    structure: str,
    output_path: str,
    *,
    node_radius: float = 0.5,
    font_size: int = 7,
    show_indices: bool = False,
    title: str | None = None,
    figsize: tuple = (10, 10),
    dpi: int = 150,
) -> str:
    sequence = sequence.upper().replace('T', 'U')
    fig, ax = plot_rna_structure(
        sequence, structure,
        node_radius=node_radius, font_size=font_size,
        show_indices=show_indices, title=title, figsize=figsize,
    )
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='none', transparent=True)
    plt.close(fig)
    return output_path