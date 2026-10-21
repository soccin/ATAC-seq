-- Table styling for md2pdf. Requires table-style.tex for \definecolor and
-- the colortbl load. LaTeX output only.
--
-- 1. Drop pandoc's computed relative column widths so LaTeX sizes each table
--    to its content instead of stretching it to \linewidth. Pandoc's markdown
--    reader assigns fractional widths only when a table's source lines are
--    wider than --columns (default 72), which makes narrow tables natural
--    width and wide ones full width in the same document. Exception: if the
--    separator dashes are unequal (|-----|------------|), keep the widths so
--    long cells wrap. Equal dashes (|---|---|) still give natural width.
--
-- 2. Shade alternate body rows, starting with the first. Done per row with
--    colortbl's \rowcolor rather than xcolor's \rowcolors, because the
--    template loads xcolor long before --include-in-header lands, so
--    [table]{xcolor} is an option clash.
--
-- 3. Drop the @{} guards from the column spec. colortbl gives every shaded
--    cell a \tabcolsep overhang on each side, so with @{} suppressing the
--    outer padding the shading is 2\tabcolsep wider than the booktabs rules.
--    Without @{} the rules gain that same padding and the two align exactly.
--    Setting the overhang to 0 instead is not equivalent: it leaves unshaded
--    gaps between the columns. Returning raw LaTeX hides the tables from
--    the writer, so md2pdf passes -V tables=true to keep the template
--    loading longtable and booktabs.

local function shade(row)
  local raw = pandoc.RawInline('latex', '\\rowcolor{tablezebra}')
  local cell = row.cells[1]
  if not cell then return end
  local blocks = cell.contents
  if #blocks > 0 and (blocks[1].t == 'Plain' or blocks[1].t == 'Para') then
    blocks[1].content:insert(1, raw)
  else
    blocks:insert(1, pandoc.Plain({raw}))
  end
end

-- True if pandoc assigned widths and they are not all equal, i.e. the author
-- set them on purpose with unequal separator dashes (|-----|------------|).
local function has_set_widths(tbl)
  local first = tbl.colspecs[1] and tbl.colspecs[1][2]
  if not first then return false end
  local unequal = false
  for _, spec in ipairs(tbl.colspecs) do
    if not spec[2] then return false end
    if math.abs(spec[2] - first) > 1e-6 then unequal = true end
  end
  return unequal
end

function Table(tbl)
  if not FORMAT:match('latex') then
    return nil
  end
  local keep_widths = has_set_widths(tbl)
  if not keep_widths then
    for _, spec in ipairs(tbl.colspecs) do
      spec[2] = nil
    end
  end
  for _, body in ipairs(tbl.bodies) do
    for i, row in ipairs(body.body) do
      if i % 2 == 1 then
        shade(row)
      end
    end
  end
  local tex = pandoc.write(pandoc.Pandoc({tbl}), 'latex')
  tex = tex:gsub('(\\begin{longtable}%[%]{)@{}(%a+)@{}}', '%1%2}')
  if keep_widths then
    -- Same @{} removal for p{} columns. Pandoc's widths subtract only the
    -- inner padding, 2(n-1)\tabcolsep; with the outer padding restored the
    -- table must subtract 2n\tabcolsep to stay within \columnwidth.
    local n = #tbl.colspecs
    tex = tex:gsub('(\\begin{longtable}%[%]{)@{}', '%1', 1)
    tex = tex:gsub('(\\real{[%d.]+}})@{}}', '%1}', 1)
    tex = tex:gsub('%(\\columnwidth %- %d*\\tabcolsep%)',
      '(\\columnwidth - ' .. 2 * n .. '\\tabcolsep)')
  end
  return pandoc.RawBlock('latex', tex)
end
