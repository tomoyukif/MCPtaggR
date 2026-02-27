#' MCPtaggR: A package to conduct mappable-collinear-polymorphic tags genotyping.
#'
#' Mappable-collinear-polymorphic tags genotyping (MCPtagg) is a pipeline for 
#' next generation sequencing (NGS)-based genotyping technique to precisely 
#' detect SNPs and count reads on validated tags (short genome segments or 
#' regions) based on genome comparison information. Short NGS reads can not 
#' always be mapped to the proper positions on the reference genome if the reads 
#' were originally derived from a genome that have plenty polymorphisms and 
#' structural variations against the reference genome, e.g. mapping reads of a 
#' hybrid population that were derived from a cross between distant relatives. 
#' Polymorphic reads that were derived from a non-reference genome could result 
#' in mismapping and unmapped due to mappings to multiple locations on the 
#' reference genome with the same mapping quality. Those mapping errors and 
#' mappability biases lead biased detection of mapped reads on SNP markers where 
#' reads of either alleles are preferentially mapped. To eliminate unexpected 
#' detection of such error prone markers, MCPtagg first filter out potential SNP 
#' markers based on the genome comparison of the founder varieties (individuals) 
#' by MUMmer followed by mappability validation of potential reads that would be 
#' sequenced by NGS. The pipeline then aligns the given reads and counts only 
#' those which were uniquely mapped onto retained mappable tags that were short 
#' genome segments simulated by in-silico digestion of the given genomes as the 
#' potential reads and then verified to have polymorphisms located in collinear 
#' blocks between the given genomes. The MCPtagg pipeline efficiently eliminates
#' error prone markers from your resultant genotype data obtained via a 
#' NGS-based genotyping technique, comapred to the result produced by the 
#' conventional pipeline that aligns reads on one reference genome.
#'
#'
#' @docType package
#' @name MCPtaggR
#' @keywords internal
"_PACKAGE"
# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
