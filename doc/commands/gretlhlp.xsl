<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for formatted GUI "online" help -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="utf-8"/>

<xsl:variable name="intl"
  select="document('hlp_l10n.xml')/internationalization"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="concat(normalize-space($intl/localization[@language=$lang]/gentext[@key=$key]/@text),' ')"/>
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

<xsl:template name="gettext-nospace">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="normalize-space($intl/localization[@language=$lang]/gentext[@key=$key]/@text)"/>
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

<xsl:template match="commandlist">
  <xsl:text>headings </xsl:text>
  <xsl:value-of select="count(command[not(@section = preceding-sibling::command/@section)])"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:for-each select="command[not(@section = preceding-sibling::command/@section)]">
    <xsl:variable name="cmd-count">
      <xsl:value-of select="count(/commandlist/command[@section = current()/@section 
                            and (not(@context) or @context=$hlp)])"/>
    </xsl:variable>
    <xsl:if test="$cmd-count &gt; 0">
      <xsl:value-of select="@section"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="$cmd-count"/>
      <xsl:text>&#10;</xsl:text>
      <xsl:for-each select="/commandlist/command[@section = current()/@section and (not(@context) or @context=$hlp)]">
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
    <xsl:text>&lt;@ref="</xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text>"&gt;&#10;</xsl:text>
  </xsl:if>
</xsl:if>
</xsl:template>

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
  <xsl:text>&#x9;</xsl:text>
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
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <xsl:text>&lt;@lit="</xsl:text>
    <xsl:value-of select="@flag"/>
    <xsl:text>"&gt;</xsl:text>
  </xsl:if> 
  <xsl:text>&lt;@var="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt; </xsl:text>
  <xsl:if test="(@optional)">] </xsl:if> 
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
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="option">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;&#x9;</xsl:text>
  </xsl:if>
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="example">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;&#x9;</xsl:text>
  </xsl:if> 
  <xsl:text>&#x9;</xsl:text>
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>  
</xsl:template>

<xsl:template match="demos">
  <xsl:choose>
    <xsl:when test="position() > 1">
      <xsl:text>&#xa;&#x9;&#x9;</xsl:text>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'Seealso'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&#x9;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="demo">
  <xsl:if test="position() > 1">
    <xsl:text>, </xsl:text>
  </xsl:if>
  <xsl:text>&lt;@inp="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
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
    <xsl:text>&#xa;&#x9;</xsl:text>
  </xsl:if> 
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="examples">
  <xsl:call-template name="nl"/>
  <xsl:choose>
    <xsl:when test="count(example) > 1">
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
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="repl">
  <xsl:text>&lt;@var="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="cmd">
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="program">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="lit|func|filename|flag">
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="mathvar">
  <xsl:text>&lt;@itl="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="super">
  <xsl:text>&lt;@sup="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="book">
  <xsl:text>&lt;@itl="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>"</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="equation">
  <xsl:if test="(@status='display')">
    <xsl:text>[/PARA]&#xa;&#xa;&#xa;  </xsl:text>
  </xsl:if>
  <xsl:text>&lt;@fig="</xsl:text>
  <xsl:value-of select="@graphic"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:if test="(@status='display')">
    <xsl:text>&#xa;&#xa;&#xa;[PARA]</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="para">
  <xsl:if test="not(@context) or @context=$hlp">
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

<xsl:template match="code">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>[CODE]</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>[/CODE]&#xa;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="pre">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp">
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

<xsl:template match="cmdref">
  <xsl:text>&lt;@ref="</xsl:text>
  <xsl:value-of select="@targ"/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="guideref">
  <xsl:text>&lt;@pdf="</xsl:text>
  <xsl:call-template name="gettext-nospace">
    <xsl:with-param name="key" select="'guidebook'"/>
  </xsl:call-template>
  <xsl:text>"&gt;</xsl:text>
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
<!-- can't handle tables at present -->
</xsl:template>

<xsl:template match="footnote"/>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

<xsl:template name="dnl">
  <xsl:text>&#10;&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
