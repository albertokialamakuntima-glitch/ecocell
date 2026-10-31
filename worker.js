// EcoCell · «Eco IA» — Cloudflare Worker que liga a Eco ao Gemini (Google).
// A chave do Gemini fica guardada aqui como segredo (GEMINI_API_KEY) e nunca aparece no site.
// O site envia a pergunta + os dados atuais (contexto); o Gemini responde em JSON: {resposta, acoes}.
// Variáveis opcionais do Worker: GEMINI_MODEL (ex.: gemini-3.6-flash) e ORIGENS_EXTRA (outros endereços autorizados, separados por vírgulas).

const ORIGENS = [
  'https://albertokialamakuntima-glitch.github.io',
  'http://localhost:8765',
  'http://127.0.0.1:8765',
];

// se um modelo não existir (404), estiver ocupado (503/500) ou sem quota (429), tenta o seguinte
const MODELOS = m => [...new Set([m, 'gemini-3.6-flash', 'gemini-3.5-flash-lite', 'gemini-flash-latest', 'gemini-flash-lite-latest', 'gemini-2.5-flash'].filter(Boolean))];

const SISTEMA = `És a Eco, a assistente da EcoCell Earth Intelligence, uma ferramenta feita para o NASA Space Apps Challenge 2026 (desafio avançado de Ciências da Terra): escolher variáveis medidas por missões ou produzidas por modelos da NASA, ver como mudam ao longo do tempo e dizer o que muda, onde, quanto e se a mudança é significativa.

Regras:
- A tua resposta vai ser LIDA EM VOZ ALTA. Escreve 1 a 4 frases curtas e naturais, sem markdown, sem listas, sem emojis, sem símbolos estranhos. Escreve os números em algarismos e arredondados (por exemplo 55, 0,2, 13 %; em inglês 0.2), nunca por extenso, e as unidades por extenso (por exemplo «graus Celsius por década»).
- Para números e factos sobre os dados usa SÓ o bloco DADOS que vem com cada pergunta. Nunca inventes valores, datas, países ou fontes.
- Se pedirem algo que os dados não têm (dados de ontem ou da semana passada quando o histórico ainda não os tem, a temperatura da água de um rio, a previsão do tempo, uma variável que não existe no site), diz com clareza que esse dado não existe no site e oferece o mais próximo que existe.
- «Ontem vs hoje», «semana passada», «atualizou?»: usa o campo «historico» (retratos diários dos dados). Diz se os números mudaram ou não e explica porquê: a NASA publica estes produtos por mês ou por ano, por isso de um dia para o outro é normal ficarem iguais. Se o histórico só tiver um dia, diz que o registo diário começou nessa data.
- Comparações entre variáveis diferentes, mesmo que não tenham nada a ver: compara o sentido da tendência, se é significativa e a percentagem da área que muda; podes usar «correlacao_series_globais», mas lembra que correlação não é causa. Diz com honestidade se a comparação faz ou não sentido físico, e porquê.
- «Significativo» = tendência de Sen testada com Mann-Kendall e corrigida com FDR (q < 0,05), em blocos de 1 grau (cerca de 110 km).
- Se a pergunta não tiver nada a ver com a Terra, a NASA ou o site, responde numa frase simpática e puxa a conversa para os dados.
- Podes mexer no mapa do site com «acoes» (no máximo 3), quando isso ajudar a mostrar a resposta:
  variavel (valor = id da variável, ex.: ndvi), pais (valor = nome em inglês, como aparece nos dados), continente (valor = nome), mundo, globo3d, mapa2d, animar (valor = id da variável, opcional), parar_animacao.
  Se não for preciso mexer no mapa, devolve acoes vazio.
- Responde na língua indicada na pergunta.`;

const ESQUEMA = {
  type: 'OBJECT',
  properties: {
    resposta: { type: 'STRING' },
    acoes: {
      type: 'ARRAY',
      items: {
        type: 'OBJECT',
        properties: {
          tipo: { type: 'STRING', enum: ['variavel', 'pais', 'continente', 'mundo', 'globo3d', 'mapa2d', 'animar', 'parar_animacao'] },
          valor: { type: 'STRING' },
        },
        required: ['tipo'],
      },
    },
  },
  required: ['resposta', 'acoes'],
};

export default {
  async fetch(req, env) {
    const origem = req.headers.get('Origin') || '';
    const extra = (env.ORIGENS_EXTRA || '').split(',').map(s => s.trim()).filter(Boolean);
    const autorizada = ORIGENS.includes(origem) || extra.includes(origem);
    const cors = {
      'Access-Control-Allow-Origin': autorizada ? origem : ORIGENS[0],
      'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
      'Access-Control-Allow-Headers': 'Content-Type',
      'Vary': 'Origin',
    };
    const json = (o, s = 200) => new Response(JSON.stringify(o), { status: s, headers: { ...cors, 'Content-Type': 'application/json; charset=utf-8' } });

    if (req.method === 'OPTIONS') return new Response(null, { status: 204, headers: cors });
    if (req.method === 'GET') return json({ ok: true, servico: 'Eco IA da EcoCell', chave: !!env.GEMINI_API_KEY });
    if (req.method !== 'POST') return json({ erro: 'método não suportado' }, 405);
    if (!autorizada) return json({ erro: 'origem não autorizada' }, 403);
    if (!env.GEMINI_API_KEY) return json({ erro: 'falta o segredo GEMINI_API_KEY no Worker' }, 500);

    let b;
    try { b = await req.json(); } catch { return json({ erro: 'pedido inválido' }, 400); }
    const pergunta = String(b.pergunta || '').slice(0, 600).trim();
    const lang = b.lang === 'en' ? 'en' : 'pt';
    const contexto = String(b.contexto || '').slice(0, 80000);
    const historico = Array.isArray(b.historico) ? b.historico.slice(-6) : [];
    if (!pergunta) return json({ erro: 'pergunta vazia' }, 400);

    // conversa anterior (para perguntas como «e no Brasil?»)
    const contents = [];
    for (const h of historico) {
      if (!h || !h.q || !h.r) continue;
      contents.push({ role: 'user', parts: [{ text: String(h.q).slice(0, 600) }] });
      contents.push({ role: 'model', parts: [{ text: JSON.stringify({ resposta: String(h.r).slice(0, 1500), acoes: [] }) }] });
    }
    contents.push({
      role: 'user',
      parts: [{ text: `DADOS (JSON, gerados agora pelo site):\n${contexto}\n\nPERGUNTA (${lang === 'en' ? 'answer in English' : 'responde em português de Portugal/Angola'}): ${pergunta}` }],
    });
    const corpo = {
      systemInstruction: { parts: [{ text: SISTEMA }] },
      contents,
      generationConfig: { temperature: 0.5, maxOutputTokens: 4096, responseMimeType: 'application/json', responseSchema: ESQUEMA },
    };

    const falhas = [];
    // pensar pouco: os dados já vêm prontos no contexto, por isso a resposta sai muito mais depressa (GEMINI_THINKING muda o nível)
    const pensar = m => /2\.5/.test(m) ? { thinkingBudget: 0 } : { thinkingLevel: env.GEMINI_THINKING || 'minimal' };
    const pede = (m, comPensar = true) => fetch(`https://generativelanguage.googleapis.com/v1beta/models/${m}:generateContent`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', 'x-goog-api-key': env.GEMINI_API_KEY },
      body: JSON.stringify(comPensar ? { ...corpo, generationConfig: { ...corpo.generationConfig, thinkingConfig: pensar(m) } } : corpo),
    });
    for (const m of MODELOS(env.GEMINI_MODEL)) {
      let r = await pede(m);
      if (r.status === 400) r = await pede(m, false);   // modelo que não aceita o nível de pensamento: pede sem ele
      if (r.status === 503 || r.status === 500) { await new Promise(ok => setTimeout(ok, 700)); r = await pede(m); }   // ocupado: mais uma tentativa
      const j = await r.json().catch(() => ({}));
      if (!r.ok) {
        const msg = (j.error && j.error.message) || `Gemini respondeu ${r.status}`;
        if ([404, 429, 500, 503].includes(r.status)) { falhas.push(`${m}: ${r.status}`); continue; }   // passa ao modelo seguinte
        return json({ erro: msg, estado: r.status }, 502);
      }
      const partes = (((j.candidates || [])[0] || {}).content || {}).parts || [];
      const texto = partes.filter(p => !p.thought).map(p => p.text || '').join('');
      let out;
      try { out = JSON.parse(texto); } catch { out = { resposta: texto, acoes: [] }; }
      return json({
        resposta: String(out.resposta || '').trim(),
        acoes: Array.isArray(out.acoes) ? out.acoes.slice(0, 3) : [],
        modelo: m,
      });
    }
    return json({ erro: `nenhum modelo disponível agora (${falhas.join('; ')})` }, 503);
  },
};
