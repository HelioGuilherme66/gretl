<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:param name="hlp">cli</xsl:param>

<xsl:output method="text" encoding="utf-8"/>

<xsl:template match="commandlist">
  <xsl:if test="$hlp='gui'">
    <xsl:text>@new-style gretl.hlp&#10;</xsl:text>     
  </xsl:if>
<xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command[not(@context) or @context=$hlp]">
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="$hlp='gui'">
      <xsl:text># </xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="@section"/>
      <xsl:text> "</xsl:text>
      <xsl:value-of select="@label"/>
      <xsl:text>"&#10;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>#&#10;</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>&#10;@</xsl:text>
      <xsl:value-of select="@section"/>      
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
  <xsl:if test="(not(@context) and $hlp='gui')">
    <xsl:text>Script command: </xsl:text>
    <xsl:value-of select="@name"></xsl:value-of>
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="command[@context and @context!=$hlp]"/>

<xsl:template match="description[not(@context) or @context=$hlp]">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="description[@context and @context!=$hlp]"/>

<xsl:template match="usage">
  <xsl:if test="$hlp='cli'">
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="arguments">
  <xsl:text>&#xa;</xsl:text>
  <xsl:choose>
    <xsl:when test="count(argument) > 1">
      <xsl:text>Arguments:  </xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>Argument:   </xsl:text> 
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:if test="(@alternate)"> or </xsl:if>
  <xsl:if test="(@optional)">[ </xsl:if> 
  <xsl:if test="@flag">
    <xsl:value-of select="@flag"/>
  </xsl:if> 
  <xsl:apply-templates/>
  <xsl:text> </xsl:text>
  <xsl:if test="(@optional)">] </xsl:if> 
</xsl:template>

<xsl:template match="options">
  <xsl:call-template name="nl"/>
  <xsl:choose>
    <xsl:when test="count(option) > 1">
      <xsl:text>Options:    </xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>Option:     </xsl:text> 
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="option|example">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;            </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altforms">
  <xsl:call-template name="nl"/>
  <xsl:text>Alternate forms:  </xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altform">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;                  </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="option|example">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;            </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="examples">
  <xsl:call-template name="nl"/>
  <xsl:choose>
    <xsl:when test="count(example) > 1">
      <xsl:text>Examples:   </xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>Example:    </xsl:text> 
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="flag">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="repl">
  <xsl:if test="@quote='true'">
    <xsl:text>"</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="@quote='true'">
    <xsl:text>"</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="cmd">
  <xsl:text>"</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="program">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="lit">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="mathvar">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="book">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>"</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="filename">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="func">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="equation">
  <xsl:if test="(@status='display')">
    <xsl:text>[/PARA]&#xa;&#xa;&#xa;  </xsl:text>
  </xsl:if>
  <xsl:value-of select="@ascii"/>
  <xsl:if test="(@status='display')">
    <xsl:text>&#xa;&#xa;&#xa;[PARA]</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="para[not(@context) or @context=$hlp]">
  <xsl:choose>
    <xsl:when test="parent::li and ancestor::ilist">
      <xsl:text>&#xa;[ILISTPAR]</xsl:text>
      <xsl:apply-templates/>[/ILISTPAR]
    </xsl:when>
    <xsl:when test="parent::li and ancestor::nlist">
      <xsl:text>&#xa;[NLISTPAR]</xsl:text>
      <xsl:apply-templates/>[/NLISTPAR]
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&#xa;[PARA]</xsl:text>
      <xsl:apply-templates/>[/PARA]
    </xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template match="para[@context and @context!=$hlp]"/>

<xsl:template match="code">
  <xsl:call-template name="dnl"/>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
</xsl:template>

<xsl:template match="ilist[not(@context) or @context=$hlp]">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="ilist[@context and @context!=$hlp]"/>

<xsl:template match="nlist[not(@context) or @context=$hlp]">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="nlist[@context and @context!=$hlp]"/>

<xsl:template match="li">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="cmdref">
<xsl:text>"</xsl:text>
  <xsl:value-of select="@targ"/>
<xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="manref">
  <xsl:text>the gretl manual</xsl:text>
</xsl:template>

<xsl:template match="tabref">
  <xsl:text>The table below </xsl:text>
</xsl:template>

<xsl:template match="menu-path">
  <xsl:if test="$hlp='cli'">
Menu path:    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="other-access">
  <xsl:if test="$hlp='cli'">
Other access: <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="table[not(@context) or @context=$hlp]"> 
<xsl:text>[TABLE]&#xa;</xsl:text>
<xsl:if test="@id">
  <xsl:text>[ROW]&#xa;[CELL]&#xa;</xsl:text>
  <xsl:value-of select="@lhead"/>
  <xsl:text>[/CELL]&#xa;[CELL]&#xa;</xsl:text>
  <xsl:value-of select="@rhead"/>
  <xsl:text>[/CELL]&#xa;[/ROW]&#xa;</xsl:text>
  <xsl:text>[ROW]&#xa;[CELL]&#xa;-------</xsl:text>
  <xsl:text>[/CELL]&#xa;[CELL]&#xa;-------</xsl:text>
  <xsl:text>[/CELL]&#xa;[/ROW]&#xa;</xsl:text>
</xsl:if>
<xsl:apply-templates/>
<xsl:text>[/TABLE]&#xa;</xsl:text>
</xsl:template>

<xsl:template match="table[@context and @context!=$hlp]"/>

<xsl:template match="row">
  <xsl:text>[ROW]&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>[/ROW]&#xa;</xsl:text>
</xsl:template>

<xsl:template match="cell">
  <xsl:text>[CELL]&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>[/CELL]&#xa;</xsl:text>
</xsl:template>

<xsl:template match="footnote"/>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

<xsl:template name="dnl">
  <xsl:text>&#10;&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
