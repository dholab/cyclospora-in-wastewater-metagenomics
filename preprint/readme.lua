local stringify = pandoc.utils.stringify

local REPOSITORY = "https://github.com/dholab/cyclospora-in-wastewater-metagenomics"
local FIGURE_WIDTHS = {
  ["02-screen-wastewater-metagenomes/results/figures/cyclospora_heatmap.svg"] = "62%",
}

local function fail(message)
  error("README preprint conversion: " .. message, 0)
end

local function is_repo_only_marker(block)
  return block.tag == "RawBlock"
    and block.format:match("html")
    and block.text:match("^%s*<!%-%-%s*repo%-only%s*%-%->%s*$")
end

local function remove_repo_only_blocks(blocks)
  local kept = pandoc.Blocks({})
  local index = 1

  while index <= #blocks do
    if is_repo_only_marker(blocks[index]) then
      local following = blocks[index + 1]
      if following == nil or following.tag ~= "Para" then
        fail("a repo-only marker must be followed by one paragraph")
      end
      index = index + 2
    else
      kept:insert(blocks[index])
      index = index + 1
    end
  end

  return kept
end

local function inline_text(inlines)
  local pieces = {}

  local function append(items)
    for _, inline in ipairs(items) do
      if inline.tag == "Str" or inline.tag == "Code" then
        pieces[#pieces + 1] = inline.text
      elseif inline.tag == "Space" then
        pieces[#pieces + 1] = " "
      elseif inline.tag == "SoftBreak" or inline.tag == "LineBreak" then
        pieces[#pieces + 1] = "\n"
      elseif inline.tag == "RawInline" then
        if not inline.format:match("html") then
          fail("unsupported raw inline in the manuscript masthead: " .. inline.format)
        end
        pieces[#pieces + 1] = inline.text
      elseif inline.tag == "Link" or inline.tag == "Emph" or inline.tag == "Strong"
          or inline.tag == "Span" or inline.tag == "Superscript" then
        append(inline.content)
      else
        fail("unsupported " .. inline.tag .. " element in the manuscript masthead")
      end
    end
  end

  append(inlines)
  return table.concat(pieces)
end

local function trim(text)
  return text:match("^%s*(.-)%s*$")
end

local function parse_authors(block)
  if block == nil or block.tag ~= "Para" then
    fail("the title must be followed by an author paragraph")
  end

  local source = inline_text(block.content)
  local authors = pandoc.MetaList({})
  local matched = 0

  for name, marks in source:gmatch("([^,]-)<sup>([^<]+)</sup>") do
    name = trim(name)
    marks = trim(marks)
    if name == "" or not marks:match("^[%d,%*]+$") then
      fail("could not understand author entry: " .. name .. " / " .. marks)
    end
    authors:insert(pandoc.MetaMap({
      name = pandoc.MetaString(name),
      marks = pandoc.MetaString(marks),
    }))
    matched = matched + 1
  end

  if matched == 0 then
    fail("no authors were found in the author paragraph")
  end

  return authors
end

local function parse_affiliations(block)
  if block == nil or block.tag ~= "Para" then
    fail("the author paragraph must be followed by affiliations")
  end

  local affiliations = pandoc.MetaList({})
  local corresponding = nil
  local line_count = 0

  for line in (inline_text(block.content) .. "\n"):gmatch("(.-)\n") do
    line = trim(line)
    if line ~= "" then
      line_count = line_count + 1
      local mark, text = line:match("^<sup>(.-)</sup>%s*(.+)$")
      if mark == nil then
        fail("affiliation line has no superscript marker: " .. line)
      end

      if mark == "*" then
        local name, email = text:match("^Corresponding author:%s*(.-)%s*•%s*(%S+)$")
        if name == nil or email == nil then
          fail("could not understand the corresponding-author line")
        end
        local email_user, email_domain = email:match("^([^@]+)@([^@]+)$")
        if email_user == nil or email_domain == nil then
          fail("could not understand the corresponding-author email address")
        end
        corresponding = pandoc.MetaMap({
          name = pandoc.MetaString(name),
          email_user = pandoc.MetaString(email_user),
          email_domain = pandoc.MetaString(email_domain),
        })
      else
        if not mark:match("^%d+$") then
          fail("invalid affiliation marker: " .. mark)
        end
        affiliations:insert(pandoc.MetaMap({
          mark = pandoc.MetaString(mark),
          text = pandoc.MetaString(text),
        }))
      end
    end
  end

  if line_count == 0 or #affiliations == 0 then
    fail("no affiliations were found")
  end
  if corresponding == nil then
    fail("no corresponding author was found")
  end

  return affiliations, corresponding
end

local function strip_caption_label(paragraph, kind, expected_number)
  if paragraph == nil or paragraph.tag ~= "Para" then
    return nil
  end
  local first = paragraph.content[1]
  if first == nil or first.tag ~= "Strong" then
    return nil
  end

  local number = stringify(first.content):match("^" .. kind .. "%s+(%d+)%.")
  if number == nil then
    return nil
  end
  number = tonumber(number)
  if number ~= expected_number then
    fail(kind .. " captions must be sequential; expected " .. expected_number
      .. " but found " .. number)
  end

  local content = first.content
  if #content < 3
      or content[1].tag ~= "Str" or content[1].text ~= kind
      or content[2].tag ~= "Space"
      or content[3].tag ~= "Str" or content[3].text ~= tostring(number) .. "." then
    fail("could not remove the label from " .. kind .. " " .. number)
  end

  local strong_content = pandoc.Inlines({})
  local start = 4
  if content[start] ~= nil and content[start].tag == "Space" then
    start = start + 1
  end
  for index = start, #content do
    strong_content:insert(content[index])
  end
  if #strong_content == 0 then
    fail(kind .. " " .. number .. " has no caption title")
  end

  local caption_inlines = pandoc.Inlines({pandoc.Strong(strong_content)})
  for index = 2, #paragraph.content do
    caption_inlines:insert(paragraph.content[index])
  end

  return number, pandoc.Caption(pandoc.Blocks({pandoc.Plain(caption_inlines)}), nil)
end

local function attach_captions(blocks)
  local converted = pandoc.Blocks({})
  local index = 1
  local next_table = 1
  local next_figure = 1

  while index <= #blocks do
    local block = blocks[index]
    local table_number, table_caption = strip_caption_label(block, "Table", next_table)

    if table_number ~= nil then
      local table_block = blocks[index + 1]
      if table_block == nil or table_block.tag ~= "Table" then
        fail("Table " .. table_number .. " caption is not followed by a table")
      end
      table_block.caption = table_caption
      table_block.identifier = "table-" .. table_number
      converted:insert(table_block)
      next_table = next_table + 1
      index = index + 2
    elseif block.tag == "Figure" then
      local figure_number, figure_caption = strip_caption_label(
        blocks[index + 1], "Figure", next_figure
      )
      if figure_number == nil then
        fail("a manuscript figure is not followed by a numbered caption")
      end
      block.caption = figure_caption
      block.identifier = "figure-" .. figure_number
      converted:insert(block)
      next_figure = next_figure + 1
      index = index + 2
    elseif block.tag == "Table" then
      if not block.classes:includes("data-availability") then
        fail("a manuscript table is missing its numbered caption")
      end
      converted:insert(block)
      index = index + 1
    else
      local figure_number = strip_caption_label(block, "Figure", next_figure)
      if figure_number ~= nil then
        fail("Figure " .. figure_number .. " caption is not preceded by a figure")
      end
      converted:insert(block)
      index = index + 1
    end
  end

  return converted
end

local function remove_section(blocks, title)
  local converted = pandoc.Blocks({})
  local index = 1
  local found = false

  while index <= #blocks do
    local block = blocks[index]
    if block.tag == "Header" and stringify(block.content) == title then
      if found then
        fail("section appears more than once: " .. title)
      end
      found = true
      local level = block.level
      index = index + 1
      while index <= #blocks do
        local candidate = blocks[index]
        if candidate.tag == "Header" and candidate.level <= level then
          break
        end
        index = index + 1
      end
    else
      converted:insert(block)
      index = index + 1
    end
  end

  if not found then
    fail("section to omit was not found: " .. title)
  end
  return converted
end

local function prepare_tables(blocks)
  local section_level = nil
  local data_availability_tables = 0

  for _, block in ipairs(blocks) do
    if block.tag == "Header" then
      if section_level ~= nil and block.level <= section_level then
        section_level = nil
      end
      if stringify(block.content) == "Data availability" then
        section_level = block.level
      end
    elseif block.tag == "Table" then
      for index, colspec in ipairs(block.colspecs) do
        if colspec[1] == pandoc.AlignDefault then
          block.colspecs[index] = {pandoc.AlignLeft, colspec[2]}
        end
      end
      if section_level ~= nil then
        data_availability_tables = data_availability_tables + 1
        block.classes:insert("data-availability")
      end
    end
  end

  if data_availability_tables ~= 1 then
    fail("expected one table in Data availability but found "
      .. data_availability_tables)
  end
end

local function is_external(target)
  return target:match("^[%a][%w+.-]*:") ~= nil
    or target:match("^//") ~= nil
    or target:match("^#") ~= nil
end

local function rewrite_link(link)
  if not is_external(link.target) then
    local path, fragment = link.target:match("^([^#]*)(#.*)$")
    if path == nil then
      path = link.target
      fragment = ""
    end
    local kind = path:match("/$") and "tree" or "blob"
    link.target = REPOSITORY .. "/" .. kind .. "/main/" .. path:gsub("/$", "") .. fragment
  end
  return link
end

local function prepare_image(image)
  if is_external(image.src) then
    fail("external images are not supported in the preprint: " .. image.src)
  end

  local original = image.src
  image.src = "../../" .. original
  if FIGURE_WIDTHS[original] ~= nil then
    image.attributes.width = FIGURE_WIDTHS[original]
  end
  return image
end

local function extract_document_parts(blocks)
  if #blocks < 5 or blocks[1].tag ~= "Header" or blocks[1].level ~= 1 then
    fail("README.md must begin with one level-one title")
  end

  local title = stringify(blocks[1].content)
  local abstract_index = nil
  for index = 2, #blocks do
    local block = blocks[index]
    if block.tag == "Header" and block.level == 2 and stringify(block.content) == "Abstract" then
      abstract_index = index
      break
    end
  end
  if abstract_index == nil then
    fail("README.md has no second-level Abstract section")
  end
  if abstract_index ~= 4 then
    fail("expected exactly an author paragraph and an affiliation paragraph before Abstract")
  end

  local authors = parse_authors(blocks[2])
  local affiliations, corresponding = parse_affiliations(blocks[3])

  local next_section = nil
  for index = abstract_index + 1, #blocks do
    local block = blocks[index]
    if block.tag == "Header" and block.level <= 2 then
      next_section = index
      break
    end
  end
  if next_section == nil then
    fail("the Abstract section is not followed by another section")
  end

  local abstract = pandoc.Blocks({})
  for index = abstract_index + 1, next_section - 1 do
    abstract:insert(blocks[index])
  end
  if #abstract == 0 then
    fail("the Abstract section is empty")
  end

  local body = pandoc.Blocks({})
  for index = next_section, #blocks do
    body:insert(blocks[index])
  end

  return title, authors, affiliations, corresponding, abstract, body
end

function Pandoc(document)
  local blocks = remove_repo_only_blocks(document.blocks)
  local title, authors, affiliations, corresponding, abstract, body =
    extract_document_parts(blocks)

  body = remove_section(body, "Reproducing this analysis")
  prepare_tables(body)

  for _, block in ipairs(body) do
    if block.tag == "Header" then
      block.level = block.level - 1
    end
  end

  body = attach_captions(body)
  body = body:walk({
    Link = rewrite_link,
    Image = prepare_image,
  })

  document.meta.title = pandoc.MetaString(title)
  document.meta.authors = authors
  document.meta.affiliations = affiliations
  document.meta.corresponding = corresponding
  document.meta.abstract = pandoc.MetaBlocks(abstract)
  document.blocks = body
  return document
end
