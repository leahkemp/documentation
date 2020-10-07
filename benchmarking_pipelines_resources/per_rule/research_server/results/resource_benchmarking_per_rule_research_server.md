Benchmarking pipeline resources - per rule - research server
================
Leah Kemp
9/15/2020

## Context

This document aims to plot and analyse the results of resource
benchmarking tests for our genomic pipelines on a single research server
(Wintermute). See related docs here: [resource benchmarking per rule
research
server](../resource_benchmarking_per_rule_research_server.md)

## Threading

### Find optimal threads for each rule (diminishing return in speed)

#### human\_genomics\_pipeline

![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

#### vcf\_annotation\_pipeline

![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

### Find maximum memory usage for each rule (max USS)

#### human\_genomics\_pipeline

![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

#### vcf\_annotation\_pipeline

![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](resource_benchmarking_per_rule_research_server_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->
