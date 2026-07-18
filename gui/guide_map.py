import math
from dataclasses import dataclass, field
from typing import Optional

from PySide6.QtCore import Qt, QPointF, Signal, QRectF
from PySide6.QtGui import (
    QColor, QPen, QBrush, QPolygonF, QPainter, QFont, QMouseEvent,
)
from PySide6.QtWidgets import (
    QGraphicsView, QGraphicsScene, QGraphicsItem,
    QGraphicsRectItem, QGraphicsPolygonItem,
    QGraphicsLineItem, QGraphicsTextItem,
)


COLOR_SENSE     = QColor(166, 227, 161)     
COLOR_ANTISENSE = QColor(116, 199, 236)    
COLOR_DISQ      = QColor(243, 139, 168)     
COLOR_REGION    = QColor(220, 220, 220)   
COLOR_SELECTED  = QColor(255, 200, 0)    
COLOR_HOVER     = QColor(255, 230, 120)   
COLOR_TRACK     = QColor(180, 190, 254)   
COLOR_TICK      = QColor(100, 100, 100)   
COLOR_TICK_TEXT = QColor(80, 80, 80)      
COLOR_TICK_MINOR = QColor(140, 140, 140)  

COLOR_PAM       = QColor(250, 179, 135)    
COLOR_PAM_TEXT  = QColor(40, 40, 40)     
COLOR_NT_TEXT   = QColor(20, 20, 20)     
COLOR_GENOME_NT = QColor(60, 60, 60)      
COLOR_CUT       = QColor(231, 130, 132)    



@dataclass
class GuideRecord:
    indx: int
    start: int
    end: int
    strand: str
    sequence: str = ""
    PAM: str = ""
    gc_freq: int = 0
    doench_16: float = 0.0
    disqualified: bool = False
    restriction_sites: list = field(default_factory=list)



class ClickableGuideItem(QGraphicsRectItem):
    def __init__(
        self,
        rect: QRectF,
        idx: int,
        guide: GuideRecord,
        base_color: QColor,
        click_signal: Signal,
        nt_mode: bool = False,
        parent: Optional[QGraphicsItem] = None,
    ):
        super().__init__(rect, parent)
        self._idx = idx
        self._guide = guide
        self._base_color = base_color
        self._click_signal = click_signal
        self._is_hovered = False
        self._is_selected = False
        self._nt_mode = nt_mode
        self._spacer_bg: Optional[QGraphicsRectItem] = None
        self._pam_bg: Optional[QGraphicsRectItem] = None
        self._arrow_poly: Optional[QGraphicsPolygonItem] = None

        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.PointingHandCursor)

        tip_parts = [
            f"{'Sense' if guide.strand == '+' else 'Antisense'}",
            f"Pos: {guide.start + 1}-{guide.end}",
        ]

        tip_parts.append(f'Index : {idx}')
        if guide.sequence:
            tip_parts.append(f"Seq: {guide.sequence}")
        if guide.PAM:
            tip_parts.append(f"PAM: {guide.PAM}")
        if guide.gc_freq:
            tip_parts.append(f"GC: {guide.gc_freq}%")
        if guide.doench_16 > 0:
            tip_parts.append(f"Doench: {guide.doench_16}")
        if guide.disqualified:
            tip_parts.append("DISQUALIFIED")
        if guide.restriction_sites:
            tip_parts.append(f"Restr: {', '.join(guide.restriction_sites)}")
        self.setToolTip("\n".join(tip_parts))

        self._update_brush()


    def hoverEnterEvent(self, event):
        self._is_hovered = True
        self._update_brush()
        self.setZValue(10)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        self._is_hovered = False
        self._update_brush()
        self.setZValue(0)
        super().hoverLeaveEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._click_signal.emit(self._idx)
            event.accept()
            return
        super().mousePressEvent(event)

    def set_selected(self, selected: bool):
        self._is_selected = selected
        self._update_brush()

    def _update_brush(self):
        if self._nt_mode:
            if self._is_selected:
                self.setBrush(QBrush(COLOR_SELECTED))
            elif self._is_hovered:
                self.setBrush(QBrush(COLOR_HOVER))
            else:
                self.setBrush(QBrush(self._base_color))
            self.setPen(QPen(self._base_color.darker(140), 1))
        else:
            self.setBrush(Qt.NoBrush)
            self.setPen(Qt.NoPen)
            if self._arrow_poly is not None:
                if self._is_selected:
                    self._arrow_poly.setBrush(QBrush(COLOR_SELECTED))
                elif self._is_hovered:
                    self._arrow_poly.setBrush(QBrush(COLOR_HOVER))
                else:
                    self._arrow_poly.setBrush(QBrush(self._base_color))

    @property
    def guide(self) -> GuideRecord:
        return self._guide

    @property
    def index(self) -> int:
        return self._idx


class GuideMapView(QGraphicsView):

    guide_clicked = Signal(int)

    TRACK_Y         = 60
    TRACK_H         = 10
    TRACK_H_NT      = 18        
    ARROW_H         = 16
    ARROW_W         = 10
    MARGIN_X        = 50
    TICK_H          = 6
    TICK_H_MINOR    = 3
    LABEL_GAP       = 4

    MAX_LANES       = 8
    LANE_H          = 20
    LANE_GAP        = 2         

    ZOOM_FACTOR     = 1.25
    NT_VIEW_THRESHOLD = 200      

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scene = QGraphicsScene(self)
        self.setScene(self._scene)
        self.setRenderHint(QPainter.Antialiasing)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setMinimumHeight(180)
        self.setDragMode(QGraphicsView.ScrollHandDrag)

        self._guides: list[GuideRecord] = []
        self._guide_items: list[ClickableGuideItem] = []
        self._selected_idx: int = -1
        self._region_start = 0
        self._region_end = 0
        self._mode = 'ko'
        self._del_start = 0
        self._del_end = 0
        self._target_seq: str = ""
        self._cut_distance: int = 3

        self._zoom_level = 1.0
        self._base_px_per_nt = 1.0
        self._coord_min = 0
        self._coord_max = 0
        self._coord_range = 1

        self._panning = False
        self._pan_start: Optional[QPointF] = None


    def set_guides(
        self,
        guides: list[GuideRecord],
        region_start: int,
        region_end: int,
        mode: str = 'ko',
        del_start: int = 0,
        del_end: int = 0,
        target_seq: str = "",
        cut_distance: int = 3,
    ):
        self._guides = guides
        self._region_start = region_start
        self._region_end = region_end
        self._mode = mode
        self._del_start = del_start
        self._del_end = del_end
        self._target_seq = target_seq or ""
        self._cut_distance = cut_distance
        self._selected_idx = -1
        self._zoom_level = 1.0
        self._redraw()

    def highlight_guide(self, idx: int):
        self._selected_idx = idx
        for item in self._guide_items:
            item.set_selected(item.index == idx)

    def clear(self):
        self._guides = []
        self._guide_items = []
        self._selected_idx = -1
        self._scene.clear()

    def reset_view(self):
        self._zoom_level = 1.0
        self._redraw()


    def _calc_base_metrics(self):
        if not self._guides:
            self._coord_min = self._region_start
            self._coord_max = max(self._region_end, self._region_start + 1)
        else:
            all_starts = [g.start for g in self._guides]
            all_ends   = [g.end   for g in self._guides]
            self._coord_min = min(all_starts + [self._region_start])
            self._coord_max = max(all_ends   + [self._region_end])
        self._coord_range = max(self._coord_max - self._coord_min, 1)

        view_w = max(self.width() - 20, 600)
        self._base_px_per_nt = (view_w - 2 * self.MARGIN_X) / self._coord_range

    def _effective_px_per_nt(self) -> float:
        return self._base_px_per_nt * self._zoom_level

    def _to_x(self, coord: int) -> float:
        return self.MARGIN_X + (coord - self._coord_min) * self._effective_px_per_nt()

    def _scene_width(self) -> float:
        return 2 * self.MARGIN_X + self._coord_range * self._effective_px_per_nt()

    def _visible_coord_range(self) -> float:
        vp_w = self.viewport().width()
        if vp_w <= 0:
            vp_w = self.width()
        return vp_w / self._effective_px_per_nt()

    def _is_nt_view(self) -> bool:
        return self._visible_coord_range() < self.NT_VIEW_THRESHOLD


    def _guide_pixel_span(self, guide: GuideRecord) -> tuple[float, float]:
        if self._is_nt_view():
            pam_len = len(guide.PAM) if guide.PAM else 3
            if guide.strand == '+':
                return self._to_x(guide.start), self._to_x(guide.end + pam_len)
            else:
                return self._to_x(guide.start - pam_len), self._to_x(guide.end)
        else:
            x_start = self._to_x(guide.start)
            x_end = self._to_x(guide.end)
            return x_start - self.LANE_GAP, x_end + self.LANE_GAP

    def _assign_lanes(self) -> dict[int, int]:
        sense_lanes: list[list[tuple[float, float]]] = [[] for _ in range(self.MAX_LANES)]
        anti_lanes:  list[list[tuple[float, float]]] = [[] for _ in range(self.MAX_LANES)]
        assignments: dict[int, int] = {}

        for idx, guide in enumerate(self._guides):
            x_s, x_e = self._guide_pixel_span(guide)
            lanes = sense_lanes if guide.strand == '+' else anti_lanes
            placed = False
            for lane_idx in range(self.MAX_LANES):
                conflict = False
                for (s, e) in lanes[lane_idx]:
                    if not (x_e <= s + self.LANE_GAP or x_s >= e - self.LANE_GAP):
                        conflict = True
                        break
                if not conflict:
                    lanes[lane_idx].append((x_s, x_e))
                    assignments[idx] = lane_idx
                    placed = True
                    break
            if not placed:
                assignments[idx] = self.MAX_LANES - 1

        return assignments


    def _calc_lane_usage(self, lane_assign: dict[int, int]) -> tuple[int, int]:
        max_sense = 0
        max_anti = 0
        for idx, guide in enumerate(self._guides):
            lane = lane_assign.get(idx, 0)
            if guide.strand == '+':
                max_sense = max(max_sense, lane + 1)
            else:
                max_anti  = max(max_anti,  lane + 1)
        return max_sense, max_anti

    def _calc_y_offset(self, max_sense: int) -> float:
        topmost_y = self.TRACK_Y - 4 - max_sense * self.LANE_H
        if topmost_y < 10:
            return 10 - topmost_y
        return 0.0

    def _calc_scene_height(self, max_sense: int, max_anti: int, y_offset: float) -> float:
        track_h = self.TRACK_H_NT if self._is_nt_view() else self.TRACK_H
        top    = self.TRACK_Y + y_offset - max_sense * self.LANE_H
        bottom = self.TRACK_Y + y_offset + track_h + max_anti * self.LANE_H
        tick_area = self.TICK_H + 20
        return bottom - top + tick_area


    def _nice_step(self, visible_range: float, target_ticks: int = 8) -> int:
        if visible_range <= 0:
            return 1
        raw = visible_range / target_ticks
        mag = 10 ** math.floor(math.log10(max(raw, 1)))
        norm = raw / mag
        if norm < 1.5:
            nice = 1
        elif norm < 3:
            nice = 2
        elif norm < 7:
            nice = 5
        else:
            nice = 10
        return max(int(nice * mag), 1)


    def _guide_color(self, guide: GuideRecord) -> QColor:
        if guide.disqualified:
            return COLOR_DISQ
        if guide.strand == '+':
            return COLOR_SENSE
        return COLOR_ANTISENSE


    def _lane_y(self, strand: str, lane: int, y_offset: float) -> float:
        track_h = self.TRACK_H_NT if self._is_nt_view() else self.TRACK_H
        if strand == '+':
            return self.TRACK_Y + y_offset - 4 - (lane + 1) * self.LANE_H
        else:
            return self.TRACK_Y + y_offset + track_h + 4 + lane * self.LANE_H


    def _redraw(self):
        self._scene.clear()
        self._guide_items = []

        if not self._guides:
            return

        self._calc_base_metrics()

        nt_view = self._is_nt_view()
        track_h = self.TRACK_H_NT if nt_view else self.TRACK_H

        lane_assign = self._assign_lanes()
        max_sense, max_anti = self._calc_lane_usage(lane_assign)
        y_offset = self._calc_y_offset(max_sense)
        self._y_offset = y_offset
        scene_h = self._calc_scene_height(max_sense, max_anti, y_offset)
        scene_w = self._scene_width()

        self._scene.setSceneRect(0, 0, scene_w, scene_h)

        track_y = self.TRACK_Y + y_offset


        track = QGraphicsRectItem(
            self.MARGIN_X, track_y,
            scene_w - 2 * self.MARGIN_X, track_h
        )
        track.setBrush(QBrush(COLOR_TRACK))
        track.setPen(QPen(Qt.NoPen))
        track.setZValue(-2)
        self._scene.addItem(track)

        if nt_view and self._target_seq:
            self._draw_genome_letters(track_h, track_y)


        if self._mode == 'deletion' and self._del_start > 0 and self._del_end > self._del_start:
            x1 = self._to_x(self._del_start - 1)
            x2 = self._to_x(self._del_end)
            del_rect = QGraphicsRectItem(
                x1, track_y - 5,
                x2 - x1, track_h + 10
            )
            del_rect.setBrush(QBrush(COLOR_REGION))
            del_rect.setPen(QPen(QColor(150, 150, 150), 1))
            del_rect.setZValue(-1)
            self._scene.addItem(del_rect)

            del_label = self._scene.addText(
                f"del {self._del_end - self._del_start + 1} bp"
            )
            del_label.setDefaultTextColor(QColor(100, 100, 100))
            del_label.setFont(QFont('monospace', 7))
            label_x = (x1 + x2) / 2 - del_label.boundingRect().width() / 2
            del_label.setPos(label_x, track_y - 22)


        self._draw_ticks(track_h, scene_w, track_y)


        for idx, guide in enumerate(self._guides):
            lane = lane_assign.get(idx, 0)
            y = self._lane_y(guide.strand, lane, y_offset)
            base_color = self._guide_color(guide)
            guide_idx = getattr(guide, 'indx', idx)

            if nt_view:
                item = self._make_nt_block(guide_idx, guide, base_color, y)
            else:
                item = self._make_arrow_item(guide_idx, guide, base_color, y)

            self._guide_items.append(item)
            self._scene.addItem(item)

        if 0 <= self._selected_idx < len(self._guide_items):
            self._guide_items[self._selected_idx].set_selected(True)

        self.setSceneRect(self._scene.itemsBoundingRect().adjusted(-10, -10, 10, 10))


    def _draw_genome_letters(self, track_h: float, track_y: float):
        px_per_nt = self._effective_px_per_nt()
        font_size = max(6, int(px_per_nt * 0.55))
        font = QFont('monospace', font_size)
        font.setStyleStrategy(QFont.StyleStrategy.NoAntialias)

        coord_left = self._coord_min
        coord_right = self._coord_max
        seq_len = len(self._target_seq)

        for coord in range(coord_left, coord_right + 1):
            if coord < 0 or coord >= seq_len:
                continue
            letter = self._target_seq[coord]
            x = self._to_x(coord)
            txt = QGraphicsTextItem(letter)
            txt.setFont(font)
            txt.setDefaultTextColor(COLOR_GENOME_NT)
            txt_w = txt.boundingRect().width()
            txt_h = txt.boundingRect().height()
            txt.setPos(x + px_per_nt / 2 - txt_w / 2,track_y + track_h / 2 - txt_h / 2)
            txt.setZValue(1)
            self._scene.addItem(txt)


    def _draw_ticks(self, track_h: float, scene_w: float, track_y: float):
        px_per_nt = self._effective_px_per_nt()
        visible_range = self._visible_coord_range()
        step = self._nice_step(visible_range)

        coord_left = self._coord_min
        coord_right = self._coord_max

        tick_y = track_y + track_h
        font = QFont('monospace', 7)

        first = math.ceil(coord_left / step) * step
        coord = first
        while coord <= coord_right:
            x = self._to_x(coord)
            if x < self.MARGIN_X - 5 or x > scene_w - self.MARGIN_X + 5:
                coord += step
                continue

            tick_line = QGraphicsLineItem(
                x, tick_y,
                x, tick_y + self.TICK_H
            )
            tick_line.setPen(QPen(COLOR_TICK, 1))
            self._scene.addItem(tick_line)

            txt = QGraphicsTextItem(str(coord))
            txt.setFont(font)
            txt.setDefaultTextColor(COLOR_TICK_TEXT)
            txt_w = txt.boundingRect().width()
            txt.setPos(x - txt_w / 2, tick_y + self.TICK_H + 2)
            self._scene.addItem(txt)

            coord += step

        minor_step = max(step // 2, 1)
        if minor_step < step:
            first_minor = math.ceil(coord_left / minor_step) * minor_step
            coord = first_minor
            while coord <= coord_right:
                if coord % step != 0:
                    x = self._to_x(coord)
                    if x >= self.MARGIN_X - 5 and x <= scene_w - self.MARGIN_X + 5:
                        tick_line = QGraphicsLineItem(
                            x, tick_y,
                            x, tick_y + self.TICK_H_MINOR
                        )
                        tick_line.setPen(QPen(COLOR_TICK_MINOR, 1))
                        self._scene.addItem(tick_line)
                coord += minor_step


    def _make_arrow_item(
        self,
        idx: int,
        guide: GuideRecord,
        color: QColor,
        y: float,
    ) -> ClickableGuideItem:
        x_start = self._to_x(guide.start)
        x_end = self._to_x(guide.end)
        h = self.ARROW_H

        span = x_end - x_start
        arrowhead_w = min(max(span * 0.15, 4), self.ARROW_W)

        if guide.strand == '+':
            pts = [
                QPointF(x_start, y),
                QPointF(x_end, y),
                QPointF(x_end + arrowhead_w, y + h / 2),
                QPointF(x_end, y + h),
                QPointF(x_start, y + h),
            ]
            rect = QRectF(x_start, y, (x_end + arrowhead_w) - x_start, h + 2)
        else:
            pts = [
                QPointF(x_end, y),
                QPointF(x_end, y + h),
                QPointF(x_start, y + h),
                QPointF(x_start - arrowhead_w, y + h / 2),
                QPointF(x_start, y),
            ]
            rect = QRectF(x_start - arrowhead_w, y, (x_end + arrowhead_w) - x_start, h + 2)

        min_click_w = 10
        if rect.width() < min_click_w:
            cx = (x_start + x_end) / 2
            rect = QRectF(cx - min_click_w / 2, y, min_click_w, h + 2)

        item = ClickableGuideItem(rect, idx, guide, color, self.guide_clicked, nt_mode=False)

        poly = QPolygonF(pts)
        arrow_poly = QGraphicsPolygonItem(poly, item)
        arrow_poly.setPen(QPen(color.darker(130), 1))
        arrow_poly.setBrush(QBrush(color))
        item._arrow_poly = arrow_poly

        item._update_brush()
        return item

    def _make_nt_block(
        self,
        idx: int,
        guide: GuideRecord,
        color: QColor,
        y: float,
    ) -> ClickableGuideItem:
        px_per_nt = self._effective_px_per_nt()
        font_size = max(6, int(px_per_nt * 0.55))
        font = QFont('monospace', font_size)
        font.setStyleStrategy(QFont.StyleStrategy.NoAntialias)

        block_h = self.LANE_H - 4

        spacer_start = guide.start
        spacer_end = guide.end

        if guide.strand == '+':
            pam_start = guide.end
            pam_end = guide.end + len(guide.PAM) if guide.PAM else guide.end + 3
        else:
            pam_end = guide.start
            pam_start = guide.start - len(guide.PAM) if guide.PAM else guide.start - 3

        span_left = min(spacer_start, pam_start)
        span_right = max(spacer_end, pam_end)
        x_left = self._to_x(span_left)
        x_right = self._to_x(span_right)
        rect = QRectF(x_left, y, x_right - x_left, block_h)

        item = ClickableGuideItem(rect, idx, guide, color, self.guide_clicked, nt_mode=True)

        pam_x = self._to_x(pam_start)
        pam_w = (pam_end - pam_start) * px_per_nt
        pam_bg = QGraphicsRectItem(QRectF(pam_x, y, pam_w, block_h), item)
        pam_bg.setBrush(QBrush(COLOR_PAM))
        pam_bg.setPen(QPen(COLOR_PAM.darker(140), 1))
        pam_bg.setZValue(1)
        item._pam_bg = pam_bg

        spacer_seq = guide.sequence[:20] if guide.sequence else ""
        for i, letter in enumerate(spacer_seq):
            if guide.strand == '+':
                coord = spacer_start + i
            else:
                coord = spacer_end - 1 - i
            x = self._to_x(coord)
            txt = QGraphicsTextItem(letter)
            txt.setParentItem(item)
            txt.setFont(font)
            txt.setDefaultTextColor(COLOR_NT_TEXT)
            txt_w = txt.boundingRect().width()
            txt_h = txt.boundingRect().height()
            txt.setPos(x + px_per_nt / 2 - txt_w / 2, y + block_h / 2 - txt_h / 2)
            txt.setZValue(2)

        pam_seq = guide.PAM[:3] if guide.PAM else ""
        for i, letter in enumerate(pam_seq):
            if guide.strand == '+':
                coord = pam_start + i
            else:
                coord = pam_end - 1 - i
            x = self._to_x(coord)
            txt = QGraphicsTextItem(letter)
            txt.setParentItem(item)
            txt.setFont(font)
            txt.setDefaultTextColor(COLOR_PAM_TEXT)
            txt_w = txt.boundingRect().width()
            txt_h = txt.boundingRect().height()
            txt.setPos(x + px_per_nt / 2 - txt_w / 2, y + block_h / 2 - txt_h / 2)
            txt.setZValue(3)

        cut_dist = getattr(self, '_cut_distance', 3)
        if guide.strand == '+':
            cut_coord = spacer_end - cut_dist - 1
        else:
            cut_coord = spacer_start + cut_dist - 1
        cut_x = self._to_x(cut_coord) + px_per_nt  
        cut_w = max(2, px_per_nt * 0.15)          
        cut_rect = QGraphicsRectItem(
            QRectF(cut_x, y - 2, cut_w, block_h + 4), item
        )
        cut_rect.setBrush(QBrush(COLOR_CUT))
        cut_rect.setPen(QPen(COLOR_CUT.darker(140), 1))
        cut_rect.setZValue(4)
        item._cut_rect = cut_rect

        arrow_size = max(4, int(px_per_nt * 0.3))
        if guide.strand == '+':
            ax = self._to_x(pam_end)
            pts = [
                QPointF(ax, y),
                QPointF(ax + arrow_size, y + block_h / 2),
                QPointF(ax, y + block_h),
            ]
        else:
            ax = self._to_x(pam_start)
            pts = [
                QPointF(ax, y),
                QPointF(ax - arrow_size, y + block_h / 2),
                QPointF(ax, y + block_h),
            ]
        arrow_poly = QGraphicsPolygonItem(QPolygonF(pts), item)
        arrow_poly.setBrush(QBrush(color.darker(130)))
        arrow_poly.setPen(QPen(color.darker(160), 1))
        arrow_poly.setZValue(2)
        item._arrow_poly = arrow_poly

        item._update_brush()
        return item


    def wheelEvent(self, event):
        if event.angleDelta().y() > 0:
            factor = self.ZOOM_FACTOR
        else:
            factor = 1 / self.ZOOM_FACTOR

        cursor_scene_pos = self.mapToScene(event.position().toPoint())
        cursor_coord = self._coord_min + (cursor_scene_pos.x() - self.MARGIN_X) / self._effective_px_per_nt()

        new_zoom = self._zoom_level * factor
        max_zoom = 30.0 / self._base_px_per_nt
        min_zoom = 0.5
        new_zoom = max(min_zoom, min(new_zoom, max_zoom))
        self._zoom_level = new_zoom

        self._redraw()

        new_x = self._to_x(cursor_coord)
        new_view_x = new_x - event.position().x()
        self.horizontalScrollBar().setValue(int(new_view_x))


    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.MiddleButton or (
            event.button() == Qt.LeftButton and event.modifiers() & Qt.ControlModifier
        ):
            self._panning = True
            self._pan_start = event.position()
            self.setCursor(Qt.ClosedHandCursor)
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent):
        if self._panning and self._pan_start is not None:
            delta = event.position() - self._pan_start
            self._pan_start = event.position()
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value() - int(delta.x())
            )
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent):
        if self._panning:
            self._panning = False
            self._pan_start = None
            self.setCursor(Qt.ArrowCursor)
            event.accept()
        else:
            super().mouseReleaseEvent(event)


    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._guides:
            self._redraw()
