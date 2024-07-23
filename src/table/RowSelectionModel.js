import {Slick} from "./Slick.js";
    // register namespace
    function RowSelectionModel(options) {
      let _grid;
      let _ranges = [];
      const _self = this;
      const _handler = new Slick.EventHandler();
      let _inHandler;
      let _options;
      const _defaults = {
        selectActiveRow: true
      };
  
      function init(grid) {
        _options = Object.assign(_defaults, options);
        _grid = grid;
        _handler.subscribe(_grid.onActiveCellChanged,
            wrapHandler(handleActiveCellChange));
        _handler.subscribe(_grid.onKeyDown,
            wrapHandler(handleKeyDown));
        _handler.subscribe(_grid.onClick,
            wrapHandler(handleClick));
      }
  
      function destroy() {
        _handler.unsubscribeAll();
      }
  
      function wrapHandler(handler) {
        return function () {
          if (!_inHandler) {
            _inHandler = true;
            handler.apply(this, arguments);
            _inHandler = false;
          }
        };
      }
  
      function rangesToRows(ranges) {
        const rows = [];
        for (let i = 0; i < ranges.length; i++) {
          for (let j = ranges[i].fromRow; j <= ranges[i].toRow; j++) {
            rows.push(j);
          }
        }
        return rows;
      }
  
      function rowsToRanges(rows) {
        const ranges = [];
        const lastCell = _grid.getColumns().length - 1;
        for (let i = 0; i < rows.length; i++) {
          ranges.push(new Slick.Range(rows[i], 0, rows[i], lastCell));
        }
        return ranges;
      }
  
      function getRowsRange(from, to) {
        let i;
        const rows = [];
        for (i = from; i <= to; i++) {
          rows.push(i);
        }
        for (i = to; i < from; i++) {
          rows.push(i);
        }
        return rows;
      }
  
      function getSelectedRows() {
        return rangesToRows(_ranges);
      }
  
      function setSelectedRows(rows) {
        setSelectedRanges(rowsToRanges(rows));
      }
  
      function setSelectedRanges(ranges) {
        // simple check for: empty selection didn't change, prevent firing onSelectedRangesChanged
        if ((!_ranges || _ranges.length === 0) && (!ranges || ranges.length === 0)) { return; }
        _ranges = ranges;
        _self.onSelectedRangesChanged.notify(_ranges);
      }
  
      function getSelectedRanges() {
        return _ranges;
      }
  
      function handleActiveCellChange(e, data) {
        if (_options.selectActiveRow && data.row != null) {
          setSelectedRanges([new Slick.Range(data.row, 0, data.row, _grid.getColumns().length - 1)]);
        }
      }
  
      function handleKeyDown(e) {
        const activeRow = _grid.getActiveCell();
        if (_grid.getOptions().multiSelect && activeRow 
        && e.shiftKey && !e.ctrlKey && !e.altKey && !e.metaKey 
        && (e.which === Slick.keyCode.UP || e.which === Slick.keyCode.DOWN)) {
          let selectedRows = getSelectedRows();
          selectedRows.sort((x, y) => x - y);
  
          if (!selectedRows.length) {
            selectedRows = [activeRow.row];
          }
  
          let top = selectedRows[0];
          let bottom = selectedRows[selectedRows.length - 1];
          let active;
  
          if (e.which === Slick.keyCode.DOWN) {
            active = activeRow.row < bottom || top === bottom ? ++bottom : ++top;
          } else {
            active = activeRow.row < bottom ? --bottom : --top;
          }
  
          if (active >= 0 && active < _grid.getDataLength()) {
            _grid.scrollRowIntoView(active);
            const tempRanges = rowsToRanges(getRowsRange(top, bottom));
            setSelectedRanges(tempRanges);
          }
  
          e.preventDefault();
          e.stopPropagation();
        }
      }
  
      function handleClick(e) {
        const cell = _grid.getCellFromEvent(e);
        if (!cell || !_grid.canCellBeActive(cell.row, cell.cell)) {
          return false;
        }
  
        if (!_grid.getOptions().multiSelect || (
            !e.ctrlKey && !e.shiftKey && !e.metaKey)) {
          return false;
        }
  
        let selection = rangesToRows(_ranges);
        const idx =  selection.indexOf(cell.row);
  
        if (idx === -1 && (e.ctrlKey || e.metaKey)) {
          selection.push(cell.row);
          _grid.setActiveCell(cell.row, cell.cell);
        } else if (idx !== -1 && (e.ctrlKey || e.metaKey)) {
          selection = selection.filter((o, i) => (o !== cell.row));
          _grid.setActiveCell(cell.row, cell.cell);
        } else if (selection.length && e.shiftKey) {
          const last = selection.pop();
          const from = Math.min(cell.row, last);
          const to = Math.max(cell.row, last);
          selection = [];
          for (let i = from; i <= to; i++) {
            if (i !== last) {
              selection.push(i);
            }
          }
          selection.push(last);
          _grid.setActiveCell(cell.row, cell.cell);
        }
  
        const tempRanges = rowsToRanges(selection);
        setSelectedRanges(tempRanges);
        e.stopImmediatePropagation();
  
        return true;
      }
  
      Object.assign(this, {
        "getSelectedRows": getSelectedRows,
        "setSelectedRows": setSelectedRows,
  
        "getSelectedRanges": getSelectedRanges,
        "setSelectedRanges": setSelectedRanges,
  
        "init": init,
        "destroy": destroy,
        "pluginName": "RowSelectionModel",
  
        "onSelectedRangesChanged": new Slick.Event()
      });
    }

    export {RowSelectionModel};
 