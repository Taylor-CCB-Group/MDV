const test_config = {
  assemblies: [
    {
      name: 'hg38',
      sequence: {
        type: 'ReferenceSequenceTrack',
        trackId: 'GRCh38-ReferenceSequenceTrack',
        adapter: {
          type: 'BgzipFastaAdapter',
          uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz',
        },
      },
      refNameAliases: {
        adapter: {
          type: 'RefNameAliasAdapter',
          uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt',
        },
      },
    },
  ],
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'genes',
      name: 'NCBI RefSeq Genes',
      assemblyNames: ['hg38'],
      category: ['Genes'],
      adapter: {
        type: 'Gff3TabixAdapter',
        uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/ncbi_refseq/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff.gz',
      },
      textSearching: {
        textSearchAdapter: {
          type: 'TrixTextSearchAdapter',
          textSearchAdapterId: 'gff3tabix_genes-index',
          uri: 'https://jbrowse.org/genomes/GRCh38/ncbi_refseq/trix/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff.gz.ix',
          assemblyNames: ['GRCh38'],
        },
      },
    },
    {
      type: 'FeatureTrack',
      trackId: 'repeats_hg38',
      name: 'Repeats',
      assemblyNames: ['hg38'],
      category: ['Annotation'],
      adapter: {
        type: 'BigBedAdapter',
        uri: 'https://jbrowse.org/genomes/GRCh38/repeats.bb',
      },
    },
    {
      type: 'AlignmentsTrack',
      trackId: 'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome',
      name: 'NA12878 Exome',
      assemblyNames: ['hg38'],
      category: ['1000 Genomes', 'Alignments'],
      adapter: {
        type: 'CramAdapter',
        uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/alignments/NA12878/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram',

        sequenceAdapter: {
          type: 'BgzipFastaAdapter',
          uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz',
        },
      },
    },
    {
      type: 'VariantTrack',
      trackId:
        'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf',
      name: '1000 Genomes Variant Calls',
      assemblyNames: ['hg38'],
      category: ['1000 Genomes', 'Variants'],
      adapter: {
        type: 'VcfTabixAdapter',
        uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/variants/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz',
      },
    },
    {
      type: 'QuantitativeTrack',
      trackId: 'hg38.100way.phyloP100way',
      name: 'hg38.100way.phyloP100way',
      category: ['Conservation'],
      assemblyNames: ['hg38'],
      adapter: {
        type: 'BigWigAdapter',
        uri: 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw',
      },
    },
    {
      type: 'AlignmentsTrack',
      trackId: 'skbr3_pacbio',
      name: 'SKBR3 pacbio',
      assemblyNames: ['hg38'],
      adapter: {
        type: 'BamAdapter',
        uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/skbr3/SKBR3_Feb17_GRCh38.sorted.bam',
      },
    },
  ],
  defaultSession: {
    name: 'this session',
    margin: 0,
    views: [
      {
        id: 'linearGenomeView',
        minimized: false,
        type: 'LinearGenomeView',
        init: {
          loc: '10:29,838,565..29,838,850',
          assembly: 'hg38',
          tracks: [
            'GRCh38-ReferenceSequenceTrack',
            'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome',
            'hg38.100way.phyloP100way',
            'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf',
          ],
        },
      },
    ],
  },
}



const assemblies={
    "hg38":{
        name: 'GRCh38',
        aliases: ['hg38'],
        sequence: {
            type: 'ReferenceSequenceTrack',
            trackId: 'GRCh38-ReferenceSequenceTrack',
            adapter: {
                type: 'BgzipFastaAdapter',
                fastaLocation: {
                    uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz',
                },
                faiLocation: {
                    uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz.fai',
                },
                gziLocation: {
                    uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz.gzi',
                },
            },
        },
        refNameAliases: {
            adapter: {
            type: 'RefNameAliasAdapter',
            location: {
                uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt',
            },
            },
        },
    },
}

const defaultView = {
    id: '_base_view',
    minimized: false,
    type: 'LinearGenomeView',
    offsetPx:85757,
    bpPerPx: 0.1554251851851852,
    displayedRegions: [
        {
        refName: '10',
        start: 0,
        end: 133797422,
        reversed: false,
        assemblyName: 'GRCh38',
        },
    ],
    tracks: [{
        id: "_base_track",
        type: "FeatureTrack",
        configuration: "_base_track",
      
        minimized: false, 
    
        displays: [{
            type: "LinearBasicDisplay",
            "heightPreConfig": 46,
            configuration: "_base_track_display",
        }]
    }],
    hideHeader: false,
    hideHeaderOverview: false,
    hideNoTracksActive: false,
    trackSelectorType: 'hierarchical',
    trackLabels: 'overlapping',
    showCenterLine: false,
    showCytobandsSetting: true,
    showGridlines: true,
}



const featureTrack= {
    type: 'FeatureTrack',
    formatDetails:{
        feature:"jexl:formatFeature(feature)"
    },
     
    trackId: '_base_track',
    name: 'NCBI RefSeq Genes',
    assemblyNames: ['GRCh38'],
    adapter: {
      type: 'BedTabixAdapter',
      bedGzLocation: {
        uri: 'tracks/genes.bed.gz',
      },
      index: {
        location: {
          uri: 'tracks/genes.bed.gz.tbi',
        },
      }
    },
    displays: [{
        type: "LinearBasicDisplay",
        displayId: "_base_track_display",
        height:100,
        renderer: {
            type: "SvgFeatureRenderer",
            color1: "jexl:colorFeature(feature)",
            height:"jexl:filterFeature(feature)",
            labels:{
                "name": "jexl:formatName(feature)"
             }
          }
        }
      ]
  }

function copy(obj){
    return JSON.parse(JSON.stringify(obj));
}

const defaultTrackConfigs={
    "bam":{
        type: 'AlignmentsTrack',   
        adapter: {
            type: 'BamAdapter',
            bamLocation: {
                uri: 'https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/jdalglei/pacbio-mm2-bam/hg38/m64176e_230414_121509.hifi_reads.mm2ax_aligned_sorted.bam',
                },
            index: {
                location:{
                    uri: 'https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/jdalglei/pacbio-mm2-bam/hg38/m64176e_230414_121509.hifi_reads.mm2ax_aligned_sorted.bam.bai',
                }
            },
            "fetchSizeLimit": 20000000,
        }
    }
    ,
    "cram":{
        type: 'AlignmentsTrack',
        adapter: {
            type: 'CramAdapter',
            cramLocation: {
              uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/alignments/NA12878/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram',
            },
            craiLocation: {
              uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/alignments/NA12878/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai',
            },
            "fetchSizeLimit": 40000000,
        }
    }

}

const defaultDisplays={
    "bam":{    
        id: 'FinKswChSr',
        displayId:"_base_alignment_display",
        type: 'LinearAlignmentsDisplay',
        PileupDisplay: {
            id: 'YAAaF494z',
            type: 'LinearPileupDisplay',
            height: 100,
            configuration: {
                type: 'LinearPileupDisplay',
                displayId:
                    'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay_LinearPileupDisplay_xyz',
            },
            showSoftClipping: true,
            filterBy: {
            flagInclude: 0,
            flagExclude: 1540,
            },
        },
        SNPCoverageDisplay: {
            id: 'VTQ_VGbAVJ',
            type: 'LinearSNPCoverageDisplay',
            height: 45,
            configuration: {
                type: 'LinearSNPCoverageDisplay',
                displayId:
                    'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay_snpcoverage_xyz',
            },
            selectedRendering: '',
            resolution: 1,
            constraints: {},
            filterBy: {
            flagInclude: 0,
            flagExclude: 1540,
            },
        },
        snpCovHeight: 45,
        height: 179,
        lowerPanelType: 'LinearPileupDisplay',
    },
}
defaultDisplays["cram"]=defaultDisplays["bam"];

const bamTrack={
    type: 'AlignmentsTrack',
    trackId: '_base_alignment',
    name: 'NA12878 Exome',
    assemblyNames: ['GRCh38'],
    
    adapter: {
        type: 'BamAdapter',
        bamLocation: {
            uri: 'https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/jdalglei/pacbio-mm2-bam/hg38/m64176e_230414_121509.hifi_reads.mm2ax_aligned_sorted.bam',
            },
        index: {
            location:{
                uri: 'https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/jdalglei/pacbio-mm2-bam/hg38/m64176e_230414_121509.hifi_reads.mm2ax_aligned_sorted.bam.bai',
            }
        },
        "fetchSizeLimit": 20000000,
        sequenceAdapter: {
            type: 'BgzipFastaAdapter',
            fastaLocation: {
                uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz',
            },
            faiLocation: {
                uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz.fai',
            },
            gziLocation: {
                uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz.gzi',
            },
        },
    },
    displays: [
        {
        id: 'FinKswChSr',
        displayId:"_base_alignment_display",
        type: 'LinearAlignmentsDisplay',
        PileupDisplay: {
            id: 'YAAaF494z',
            type: 'LinearPileupDisplay',
            height: 100,
            configuration: {
            type: 'LinearPileupDisplay',
            displayId:
                'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay_LinearPileupDisplay_xyz',
            },
            showSoftClipping: true,
            filterBy: {
            flagInclude: 0,
            flagExclude: 1540,
            },
        },
        SNPCoverageDisplay: {
            id: 'VTQ_VGbAVJ',
            type: 'LinearSNPCoverageDisplay',
            height: 45,
            configuration: {
            type: 'LinearSNPCoverageDisplay',
            displayId:
                'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay_snpcoverage_xyz',
            },
            selectedRendering: '',
            resolution: 1,
            constraints: {},
            filterBy: {
            flagInclude: 0,
            flagExclude: 1540,
            },
        },
        snpCovHeight: 45,
        height: 179,
        lowerPanelType: 'LinearPileupDisplay',
        },
    ],
}





const bamTrackView= {
    id: '_def_alignment',
    type: 'AlignmentsTrack',
    configuration: '_def_alignment',
    minimized: false,
    displays: [
        {
        id: 'FinKswChSr',
        type: 'LinearAlignmentsDisplay',
        PileupDisplay: {
            id: 'YAAaF494z',
            type: 'LinearPileupDisplay',
            height: 134,
            configuration: {
            type: 'LinearPileupDisplay',
            displayId:
                'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay_LinearPileupDisplay_xyz',
            },
            showSoftClipping: false,
            filterBy: {
            flagInclude: 0,
            flagExclude: 1540,
            },
        },
        SNPCoverageDisplay: {
            id: 'VTQ_VGbAVJ',
            type: 'LinearSNPCoverageDisplay',
            height: 45,
            configuration: {
            type: 'LinearSNPCoverageDisplay',
            displayId:
                'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay_snpcoverage_xyz',
            },
            selectedRendering: '',
            resolution: 1,
            constraints: {},
            filterBy: {
            flagInclude: 0,
            flagExclude: 1540,
            },
        },
        snpCovHeight: 45,
        configuration:
            'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome-LinearAlignmentsDisplay',
        height: 179,
        lowerPanelType: 'LinearPileupDisplay',
        },
    ],
}

function getJBrowseConfig(genome="hg38",defaultTrack,{BPView=false}){
    const def_track = copy(featureTrack);
    def_track.assemblyNames=[genome];
    def_track.name = defaultTrack.label;
    def_track.height=100;
    def_track.adapter.bedGzLocation.uri = defaultTrack.url;
    def_track.adapter.index.location.uri = defaultTrack.url+".tbi";
    let view = copy(defaultView);
    if (BPView){
        view={
            
            "id": "Dg4C2QrtR_N1qNt7CX-OT",
            "displayName": "breakend split detail",
            "minimized": false,
            "type": "BreakpointSplitView",
            "height": 400,
            "trackSelectorType": "hierarchical",
            "showIntraviewLinks": true,
            "linkViews": false,
            "interactToggled": false,
            views:[view]
        }
    }
    const conf= {
        assemblies:[copy(assemblies[genome])],
        tracks:[copy(def_track)],
        defaultSession:{
            name:"Default",
            margin:0,
            views:[view]
        }
    }
    return conf;
}


function getTrackConfig({genome,url,name,id,type}){
    const track = copy(defaultTrackConfigs[type]);
    track.assemblyNames=[genome];
    track.name = name
    track.trackId = id;
    if (type==="bam"){
        track.adapter.bamLocation.uri = url;
        track.adapter.index.location.uri = url+".bai";
        track.adapter.sequenceAdapter = copy(assemblies[genome].sequence.adapter);
    }
    if (type==="cram"){
        track.adapter.cramLocation.uri = url;
        track.adapter.craiLocation.uri = url+".crai";
        track.adapter.sequenceAdapter = copy(assemblies[genome].sequence.adapter);
    }

    const display = copy(defaultDisplays[type]);
    return [track,display];
}

export {getJBrowseConfig,getTrackConfig,test_config}