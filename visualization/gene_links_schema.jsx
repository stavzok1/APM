import React, { useState } from 'react';

const SchemaNode = ({ label, type, children, isKey = false, isRoot = false, depth = 0 }) => {
  const [expanded, setExpanded] = useState(depth < 3);
  const hasChildren = children && children.length > 0;
  
  const colors = {
    dict: '#3b82f6',      // blue
    list: '#8b5cf6',      // purple
    string: '#10b981',    // green
    number: '#f59e0b',    // amber
    boolean: '#ef4444',   // red
    key: '#6b7280',       // gray
    root: '#1f2937',      // dark
  };
  
  const bgColors = {
    dict: '#eff6ff',
    list: '#f5f3ff',
    string: '#ecfdf5',
    number: '#fffbeb',
    boolean: '#fef2f2',
  };

  const getTypeIcon = (t) => {
    switch(t) {
      case 'dict': return '{ }';
      case 'list': return '[ ]';
      case 'string': return 'str';
      case 'number': return '#';
      case 'boolean': return 'T/F';
      default: return '?';
    }
  };

  return (
    <div className={`${isRoot ? '' : 'ml-4'} my-1`}>
      <div 
        className={`flex items-center gap-2 p-2 rounded-lg cursor-pointer transition-all hover:shadow-md ${hasChildren ? 'hover:bg-gray-50' : ''}`}
        onClick={() => hasChildren && setExpanded(!expanded)}
        style={{ 
          borderLeft: isRoot ? 'none' : `3px solid ${colors[type] || colors.key}`,
          backgroundColor: isRoot ? '#f8fafc' : (expanded && hasChildren ? bgColors[type] : 'white')
        }}
      >
        {hasChildren && (
          <span className="text-gray-400 w-4 text-center font-mono">
            {expanded ? '▼' : '▶'}
          </span>
        )}
        {!hasChildren && <span className="w-4" />}
        
        <span 
          className="px-2 py-0.5 rounded text-xs font-mono font-bold"
          style={{ 
            backgroundColor: colors[type] || colors.key,
            color: 'white'
          }}
        >
          {getTypeIcon(type)}
        </span>
        
        <span className={`font-medium ${isKey ? 'text-gray-500 italic' : 'text-gray-800'}`}>
          {label}
        </span>
        
        {isKey && (
          <span className="text-xs text-gray-400 ml-1">(dynamic key)</span>
        )}
      </div>
      
      {expanded && hasChildren && (
        <div className="border-l-2 border-gray-200 ml-3">
          {children}
        </div>
      )}
    </div>
  );
};

const Legend = () => (
  <div className="flex flex-wrap gap-4 p-4 bg-gray-50 rounded-lg mb-6">
    <div className="text-sm font-semibold text-gray-600">Legend:</div>
    {[
      { type: 'dict', label: 'Dictionary { }' },
      { type: 'list', label: 'List [ ]' },
      { type: 'string', label: 'String' },
      { type: 'number', label: 'Number' },
      { type: 'boolean', label: 'Boolean' },
    ].map(({ type, label }) => (
      <div key={type} className="flex items-center gap-1">
        <span 
          className="px-2 py-0.5 rounded text-xs font-mono font-bold text-white"
          style={{ backgroundColor: type === 'dict' ? '#3b82f6' : type === 'list' ? '#8b5cf6' : type === 'string' ? '#10b981' : type === 'number' ? '#f59e0b' : '#ef4444' }}
        >
          {type === 'dict' ? '{ }' : type === 'list' ? '[ ]' : type === 'string' ? 'str' : type === 'number' ? '#' : 'T/F'}
        </span>
        <span className="text-sm text-gray-600">{label}</span>
      </div>
    ))}
    <div className="flex items-center gap-1">
      <span className="text-gray-500 italic text-sm">italic</span>
      <span className="text-sm text-gray-600">= dynamic key</span>
    </div>
  </div>
);

const GeneLinkSchema = () => {
  return (
    <div className="p-6 max-w-4xl mx-auto bg-white min-h-screen">
      <h1 className="text-2xl font-bold text-gray-800 mb-2">gene_links Schema Structure</h1>
      <p className="text-gray-600 mb-4">
        Nested data structure for enhancer-gene regulatory links from ENCODE SCREEN and ABC model
      </p>
      
      <Legend />
      
      <div className="border rounded-xl p-4 bg-white shadow-sm">
        <SchemaNode label="gene_links" type="dict" isRoot={true} depth={0}>
          <SchemaNode label="<gene_name>" type="dict" isKey={true} depth={1}>
            
            {/* Gene metadata */}
            <SchemaNode label="gene_id" type="string" depth={2} />
            <SchemaNode label="gene_type" type="string" depth={2} />
            <SchemaNode label="region" type="string" depth={2} />
            
            {/* SCREEN Experimental */}
            <SchemaNode label="screen_exp" type="dict" depth={2}>
              <SchemaNode label="per_biosample" type="dict" depth={3}>
                <SchemaNode label="<biosample_name>" type="dict" isKey={true} depth={4}>
                  <SchemaNode label="<assay_type>" type="dict" isKey={true} depth={5}>
                    <SchemaNode label="score" type="number" depth={6} />
                    <SchemaNode label="p_value" type="number" depth={6} />
                    <SchemaNode label="strength" type="string" depth={6} />
                  </SchemaNode>
                </SchemaNode>
              </SchemaNode>
              <SchemaNode label="conservation_global" type="dict" depth={3}>
                <SchemaNode label="<assay_type>" type="dict" isKey={true} depth={4}>
                  <SchemaNode label="n_biosamples" type="number" depth={5} />
                  <SchemaNode label="n_strong" type="number" depth={5} />
                  <SchemaNode label="n_weak" type="number" depth={5} />
                </SchemaNode>
              </SchemaNode>
              <SchemaNode label="conservation_breast" type="dict" depth={3}>
                <SchemaNode label="<assay_type>" type="dict" isKey={true} depth={4}>
                  <SchemaNode label="n_biosamples" type="number" depth={5} />
                  <SchemaNode label="n_strong" type="number" depth={5} />
                  <SchemaNode label="n_weak" type="number" depth={5} />
                </SchemaNode>
              </SchemaNode>
            </SchemaNode>
            
            {/* SCREEN Computational */}
            <SchemaNode label="screen_comp" type="dict" depth={2}>
              <SchemaNode label="per_biosample" type="dict" depth={3}>
                <SchemaNode label="<biosample_name>" type="dict" isKey={true} depth={4}>
                  <SchemaNode label="<assay_type>" type="dict" isKey={true} depth={5}>
                    <SchemaNode label="score" type="number" depth={6} />
                    <SchemaNode label="p_value" type="number" depth={6} />
                    <SchemaNode label="strength" type="string" depth={6} />
                  </SchemaNode>
                </SchemaNode>
              </SchemaNode>
              <SchemaNode label="conservation_global" type="dict" depth={3}>
                <SchemaNode label="<assay_type>" type="dict" isKey={true} depth={4}>
                  <SchemaNode label="n_biosamples" type="number" depth={5} />
                  <SchemaNode label="n_strong" type="number" depth={5} />
                  <SchemaNode label="n_weak" type="number" depth={5} />
                </SchemaNode>
              </SchemaNode>
              <SchemaNode label="conservation_breast" type="dict" depth={3}>
                <SchemaNode label="<assay_type>" type="dict" isKey={true} depth={4}>
                  <SchemaNode label="n_biosamples" type="number" depth={5} />
                  <SchemaNode label="n_strong" type="number" depth={5} />
                  <SchemaNode label="n_weak" type="number" depth={5} />
                </SchemaNode>
              </SchemaNode>
            </SchemaNode>
            
            {/* ABC Enhancers */}
            <SchemaNode label="ABC_enhancers" type="list" depth={2}>
              <SchemaNode label="<enhancer>" type="dict" isKey={true} depth={3}>
                <SchemaNode label="start" type="number" depth={4} />
                <SchemaNode label="end" type="number" depth={4} />
                <SchemaNode label="ABC_full" type="dict" depth={4}>
                  <SchemaNode label="<celltype>" type="dict" isKey={true} depth={5}>
                    <SchemaNode label="ABC_score" type="number" depth={6} />
                    <SchemaNode label="ABC_num" type="number" depth={6} />
                    <SchemaNode label="activity" type="number" depth={6} />
                    <SchemaNode label="distance" type="number" depth={6} />
                    <SchemaNode label="element_class" type="string" depth={6} />
                    <SchemaNode label="is_self_promoter" type="boolean" depth={6} />
                    <SchemaNode label="hic_pl_scaled" type="number" depth={6} />
                    <SchemaNode label="powerlaw_score" type="number" depth={6} />
                    <SchemaNode label="gene_expr" type="number" depth={6} />
                    <SchemaNode label="promoter_activity_q" type="number" depth={6} />
                    <SchemaNode label="gene_is_expressed" type="boolean" depth={6} />
                    <SchemaNode label="rank_within_gene" type="number" depth={6} />
                    <SchemaNode label="is_present" type="boolean" depth={6} />
                    <SchemaNode label="is_strong" type="boolean" depth={6} />
                  </SchemaNode>
                </SchemaNode>
              </SchemaNode>
            </SchemaNode>
            
          </SchemaNode>
        </SchemaNode>
      </div>
      
      {/* Path examples */}
      <div className="mt-6 p-4 bg-gray-50 rounded-lg">
        <h3 className="font-semibold text-gray-700 mb-3">Example Query Paths</h3>
        <div className="space-y-2 font-mono text-sm">
          <div className="p-2 bg-white rounded border">
            <span className="text-blue-600">screen_exp.per_biosample.</span>
            <span className="text-purple-600">*</span>
            <span className="text-blue-600">.</span>
            <span className="text-purple-600">*</span>
            <span className="text-blue-600">.strength</span>
            <span className="text-gray-400 ml-2">→ any biosample, any assay</span>
          </div>
          <div className="p-2 bg-white rounded border">
            <span className="text-blue-600">ABC_enhancers.</span>
            <span className="text-purple-600">[*]</span>
            <span className="text-blue-600">.ABC_full.</span>
            <span className="text-purple-600">*</span>
            <span className="text-blue-600">.ABC_score</span>
            <span className="text-gray-400 ml-2">→ all enhancers, any celltype</span>
          </div>
          <div className="p-2 bg-white rounded border">
            <span className="text-blue-600">screen_exp.conservation_breast.</span>
            <span className="text-purple-600">*</span>
            <span className="text-blue-600">.n_strong</span>
            <span className="text-gray-400 ml-2">→ any assay conservation</span>
          </div>
        </div>
      </div>
      
      {/* Data source labels */}
      <div className="mt-6 grid grid-cols-3 gap-4">
        <div className="p-3 rounded-lg bg-blue-50 border border-blue-200">
          <div className="font-semibold text-blue-800">screen_exp</div>
          <div className="text-sm text-blue-600">ENCODE SCREEN 3D-Chromatin</div>
          <div className="text-xs text-blue-500 mt-1">Intact-HiC, ChIA-PET, RNAPII-ChIAPET </div>
        </div>
        <div className="p-3 rounded-lg bg-green-50 border border-green-200">
          <div className="font-semibold text-green-800">screen_comp</div>
          <div className="text-sm text-green-600">ENCODE SCREEN Computational</div>
          <div className="text-xs text-green-500 mt-1">ABC-DNase, EPIraction, rE2G</div>
        </div>
        <div className="p-3 rounded-lg bg-purple-50 border border-purple-200">
          <div className="font-semibold text-purple-800">ABC_enhancers</div>
          <div className="text-sm text-purple-600">ABC Model Predictions</div>
          <div className="text-xs text-purple-500 mt-1">Full Hi-C based predictions</div>
        </div>
      </div>
    </div>
  );
};

export default GeneLinkSchema;
