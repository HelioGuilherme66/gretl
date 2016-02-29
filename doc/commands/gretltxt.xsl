<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for gretlcli plain text "online" help -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="utf-8"/>

<xsl:variable name="intl"
	      select="document('hlp_l10n.xml')/internationalization"/>

<xsl:variable name="docref"
	      select="document('docref.xml')/refsets"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="$intl/localization[@language=$lang]/gentext[@key=$key]/@text"/>
  <xsl:choose>
    <xsl:when test="$itext">
      <xsl:value-of select="$itext"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:message terminate="yes">
        <xsl:text>** Error: no phrase with key = '</xsl:text>
        <xsl:value-of select="$key"/>
        <xsl:text>' found for lang '</xsl:text>
        <xsl:value-of select="$lang"/>
        <xsl:text>'.</xsl:text>
      </xsl:message>
    </xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template name="getref">
  <xsl:param name="key"/>
  <xsl:variable name="chap"
    select="normalize-space($docref/refset[@id='guide-chapters']/ref[@key=$key]/@chapter)"/>
  <xsl:choose>
    <xsl:when test="$chap">
      <xsl:value-of select="$chap"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:message terminate="yes">
        <xsl:text>** Error: no phrase with key = '</xsl:text>
        <xsl:value-of select="$key"/>
        <xsl:text>' found.'</xsl:text>
      </xsl:message>
    </xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template match="commandref">
  <xsl:text>headings </xsl:text>
  <xsl:value-of select="count(command[not(@section = preceding-sibling::command/@section)])"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:for-each select="command[not(@section = preceding-sibling::command/@section)]">
    <xsl:variable name="cmd-count">
      <xsl:value-of select="count(/commandref/command[@section = current()/@section 
                            and (not(@context) or @context=$hlp)])"/>
    </xsl:variable>
    <xsl:if test="$cmd-count &gt; 0">
      <xsl:value-of select="@section"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="$cmd-count"/>
      <xsl:text>&#10;</xsl:text>
      <xsl:for-each select="/commandref/command[@section = current()/@section and (not(@context) or @context=$hlp)]">
        <xsl:value-of select="@name"/>
        <xsl:if test="$hlp='gui'">
          <xsl:text> "</xsl:text>
          <xsl:value-of select="@label"/>
          <xsl:text>"</xsl:text>
        </xsl:if>
        <xsl:text>&#10;</xsl:text>
      </xsl:for-each>      
    </xsl:if>
  </xsl:for-each>
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
  <xsl:if test="not(@context) or @context=$hlp">
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:text># </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text> </xsl:text>
  <xsl:value-of select="@section"/>
  <xsl:if test="$hlp='gui'">
    <xsl:text> "</xsl:text>
    <xsl:value-of select="@label"/>
    <xsl:text>"</xsl:text>
  </xsl:if>
  <xsl:text>&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
  <xsl:if test="(not(@context) and $hlp='gui')">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'scriptcommand'"/>
    </xsl:call-template>
    <xsl:value-of select="@name"></xsl:value-of>
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:if>
</xsl:template>

<xsl:template match="common-opt" />

<xsl:template match="description">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

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
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'args'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'arg'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="argblock">
  <xsl:if test="(@optional)">[ </xsl:if>
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="(@optional)">] </xsl:if>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@optional)">[ </xsl:if> 
  <xsl:if test="(@separated and (preceding-sibling::argument or preceding-sibling::argblock))">; </xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <xsl:value-of select="@flag"/>
  </xsl:if> 
  <xsl:apply-templates/>
  <xsl:text> </xsl:text>
  <xsl:if test="(@separated and not(preceding-sibling::argument or preceding-sibling::argblock))">; </xsl:if>
  <xsl:if test="(@optional)">] </xsl:if> 
</xsl:template>

<xsl:template match="argpunct">
   <xsl:apply-templates/>
</xsl:template>

<xsl:template match="options">
  <xsl:call-template name="nl"/>
  <xsl:choose>
    <xsl:when test="count(option) > 1">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'opts'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'opt'"/>
      </xsl:call-template>
      <xsl:if test="$hlp='gui'">
       <xsl:text> </xsl:text>
      </xsl:if>
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

<xsl:template match="optparm">
  <xsl:if test="(@optional)">[</xsl:if> 
  <xsl:text>=</xsl:text>
  <xsl:apply-templates/>
  <xsl:if test="(@optional)">]</xsl:if> 
</xsl:template>

<xsl:template match="demos">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;            </xsl:text>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'Seealso'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="demo">
  <xsl:if test="position() > 1">
    <xsl:text>, </xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altforms">
  <xsl:call-template name="nl"/>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'altforms'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altform">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;            </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="syntax">
  <xsl:call-template name="nl"/>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'syntax'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="examples">
  <xsl:if test="ancestor::command">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="count(example) > 1 or count(demos) > 0">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'examples'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'example'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:if test="ancestor::function">
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="flag">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="note">
  <xsl:apply-templates/>
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

<xsl:template match="funcref">
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="funclist">
  <xsl:text>## </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="function">
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:text># </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'output'"/>
  </xsl:call-template>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="@output"/>
  </xsl:call-template>
  <xsl:if test="@altout">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="@altout"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:text>&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
</xsl:template>

<xsl:template match="fnargs">
  <xsl:choose>
    <xsl:when test="count(fnarg) > 1">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'args'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'arg'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="fnarg">
  <xsl:if test="position() > 1">
    <xsl:text>            </xsl:text>
  </xsl:if> 
  <xsl:choose>
    <xsl:when test="@type='varargs'">
      <xsl:text>...</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text> (</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="@type"/>
  </xsl:call-template>
  <xsl:if test="(@optional)">
    <xsl:choose>
      <xsl:when test="@type='matrixref'">
        <xsl:text>,</xsl:text>
        <xsl:call-template name="gettext">
          <xsl:with-param name="key" select="'or'"/>
        </xsl:call-template>
        <xsl:text>null</xsl:text>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>, </xsl:text>
        <xsl:call-template name="gettext">
          <xsl:with-param name="key" select="'optional'"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:if>
  <xsl:if test="(@conditional)">
    <xsl:text>, </xsl:text>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'conditional'"/>
    </xsl:call-template>
  </xsl:if> 
  <xsl:text>)</xsl:text>
  <xsl:text>&#10;</xsl:text>
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
  <xsl:if test="@quote='true'">
    <xsl:text>"</xsl:text> 
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="@quote='true'">
    <xsl:text>"</xsl:text>
  </xsl:if>  
</xsl:template>

<xsl:template match="opt">
  <xsl:choose>
    <xsl:when test="substring(text(),1,2)='--'">
      <xsl:apply-templates/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>--</xsl:text> 
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="math">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="sup">
  <xsl:text>^</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="sub">
  <xsl:text>_</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="book">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="cite">
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

<xsl:template match="para">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
  <xsl:choose>
    <xsl:when test="parent::li and ancestor::ilist">
      <xsl:text>&#xa;[ILISTPAR]</xsl:text>
      <xsl:apply-templates/>[/ILISTPAR]
    </xsl:when>
    <xsl:when test="parent::li and ancestor::nlist">
      <xsl:text>&#xa;[NLISTPAR]</xsl:text>
      <xsl:if test="count(../*)=1">
        <xsl:value-of select="1 + count(../preceding-sibling::*)"/>
        <xsl:text>. </xsl:text>
      </xsl:if>
      <xsl:apply-templates/>[/NLISTPAR]
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&#xa;[PARA]</xsl:text>
      <xsl:apply-templates/>[/PARA]
    </xsl:otherwise>
  </xsl:choose> 
 </xsl:if>
</xsl:template>

<xsl:template match="refnote">
  <xsl:if test="@xref='false'">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="pre">
  <xsl:if test="not(@context) or @context=$hlp">
      <xsl:text>&#xa;[PRE]</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>[/PRE]</xsl:text>
 </xsl:if>
</xsl:template>

<xsl:template match="code|mono">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:call-template name="dnl"/>
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="nlist">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="li">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="subhead">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>[PARA]</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>[/PARA]&#xa;</xsl:text>    
  </xsl:if>
</xsl:template>

<xsl:template match="cmdref">
  <xsl:text>"</xsl:text>
  <xsl:value-of select="@targ"/>
  <xsl:text>"</xsl:text> 
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>     
</xsl:template>

<xsl:template match="fncref">
  <xsl:text>"</xsl:text>
  <xsl:value-of select="@targ"/>
  <xsl:text>"</xsl:text>      
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="guideref">
  <xsl:choose>
    <xsl:when test="$lang = 'en'">
      <xsl:text>chapter </xsl:text>
      <xsl:call-template name="getref">
	<xsl:with-param name="key" select="@targ"/>
      </xsl:call-template>
      <xsl:text> of the Gretl User's Guide</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
	<xsl:with-param name="key" select="'guidebook'"/>
      </xsl:call-template>
       <xsl:text> (</xsl:text>
      <xsl:call-template name="gettext">
	<xsl:with-param name="key" select="'chapter'"/>
      </xsl:call-template>
      <xsl:text> </xsl:text>
      <xsl:call-template name="getref">
	<xsl:with-param name="key" select="@targ"/>
      </xsl:call-template>
      <xsl:text>)</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="menu-path">
  <xsl:text>&#xa;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'menupath'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="other-access">
  <xsl:text>&#xa;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'otheraccess'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="table"> 
<xsl:if test="not(@context) or @context=$hlp">
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
</xsl:if>
</xsl:template>

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

<xsl:template match="seelist">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'Seealso'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
  <xsl:text>.</xsl:text>
</xsl:template>

<xsl:template match="by">
  <xsl:value-of select="@r"/>
  <xsl:text> x </xsl:text>
  <xsl:value-of select="@c"/>
</xsl:template>

</xsl:stylesheet>
